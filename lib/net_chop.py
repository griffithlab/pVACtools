import argparse
import sys
import requests
import csv
import tempfile
import re
import os
import pymp
from time import sleep
from Bio import SeqIO
import collections
import lib.run_utils

cycle = ['|', '/', '-', '\\']
methods = ['cterm', '20s']

class NetChop:
    def __init__(self, input_file, input_fasta, output_file, method='cterm', threshold=0.5, file_type='pVACseq', n_threads=1):
        self.input_file = input_file
        self.input_fasta = input_fasta
        self.output_file = output_file
        self.method = method
        self.threshold = float(threshold)
        self.file_type = file_type
        self.n_threads = 1 if sys.platform == "darwin" else n_threads # pymp and requests not compatible on macOS 10.13+ for n_threads > 1

    def get_mt_peptides(self):
        records = list(SeqIO.parse(self.input_fasta, "fasta"))
        if self.file_type == 'pVACseq':
            records_dict = {x.id.replace('MT.', ''): str(x.seq) for x in filter(lambda x: x.id.startswith('MT.'), records)}
        else:
            records_dict = {x.id: str(x.seq) for x in records}
        return records_dict

    def extract_flanked_epitope(self, full_peptide, epitope):
        flanking_sequence_length = 9
        ep_start = full_peptide.index(epitope)
        start = ep_start - flanking_sequence_length
        if start < 0:
            start = 0
        start_diff = ep_start - start
        end = ep_start + len(epitope) + flanking_sequence_length
        return full_peptide[start:end], start_diff

    def populate_staging_file(self, chunk, staging_file, current_buffer, seqs_start_diff, mt_records_dict, seq_id):
        for line in chunk:
            sequence_id = ('%010x'%seq_id[0])[-10:]
            staging_file.write('>'+sequence_id+'\n')
            if self.file_type == 'pVACbind' or self.file_type == 'pVACfuse':
                full_peptide = mt_records_dict[line['Mutation']]
                epitope = line['Epitope Seq']
            else:
                full_peptide = mt_records_dict[line['Index']]
                epitope = line['MT Epitope Seq']
            peptide, start_diff = self.extract_flanked_epitope(full_peptide, epitope)
            staging_file.write(peptide+'\n')
            current_buffer[sequence_id] = {k:line[k] for k in line}
            seqs_start_diff[sequence_id] = (start_diff, len(epitope))
            seq_id[0]+=1
        staging_file.seek(0)

    def call_net_chop(self, staging_file, chosen_method, cycle_ind, pmp):
        jobid_searcher = re.compile(r'<!-- jobid: [0-9a-fA-F]*? status: (queued|active)')
        result_delimiter = re.compile(r'-{20,}')
        fail_searcher = re.compile(r'(Failed run|Problematic input:)')

        with pmp.lock: # stagger calls to NetChop server
            sleep(5)

        response = requests.post(
            "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi",
            files={'SEQSUB':(staging_file.name, staging_file, 'text/plain')},
            data = {
                'configfile':'/var/www/html/services/NetChop-3.1/webface.cf',
                'SEQPASTE':'',
                'method':chosen_method,
                'thresh':'%0f'%self.threshold
            }
        )

        while jobid_searcher.search(response.content.decode()):
            for _ in range(10):
                with pmp.lock:
                    sys.stdout.write('\b'+cycle[cycle_ind[0]%4])
                    cycle_ind[0]+=1
                    sys.stdout.flush()
                sleep(1)
            response = requests.get(response.url)
        if fail_searcher.search(response.content.decode()):
            sys.stdout.write('\b\b')
            print('Failed!')
            print("NetChop encountered an error during processing")
            sys.exit(1)

        results = [item.strip() for item in result_delimiter.split(response.content.decode())]
        return results

    def write_line(self, line, writer, current_chunk_id, out_chunk_id, pmp):
        printed = False
        while not printed:
            with pmp.lock:
                if current_chunk_id == out_chunk_id[0]:
                    writer.writerow(line)
                    printed = True
            if current_chunk_id != out_chunk_id[0]: sleep(1)

    def write_results(self, results, current_buffer, seqs_start_diff, writer, current_chunk_id, out_chunk_id, pmp):
        for i in range(2, len(results), 4): #examine only the parts we want, skipping all else
            score = -1
            pos = 0
            sequence_name = False
            start_diff = False
            ep_len = False
            cleavage_scores = {}
            for line in results[i].split('\n'):
                data = [word for word in line.strip().split(' ') if len(word)]
                if not sequence_name:
                    sequence_name = data[4]
                    start_diff = seqs_start_diff[sequence_name][0]
                    ep_len = seqs_start_diff[sequence_name][1]
                currentPosition = data[0]
                isCleavage = data[2]
                if isCleavage is not 'S':
                    continue
                currentScore = float(data[3])
                cleavage_scores[currentPosition] = currentScore
            if len(cleavage_scores) == 0:
                best_cleavage_position = 'NA'
                best_cleavage_score = 'NA'
                cleavage_sites = 'NA'
            else:
                #filter out cleavage sites outside epitope and adjust positions in accordance
                epitope_cleavage_scores = [
                    (x[0] - start_diff, x[1])
                    for x in map(lambda x: (int(x[0]), x[1]), cleavage_scores.items())
                    if x[0] >= start_diff and x[0] <= start_diff + ep_len
                ]
                max_cleavage_score = max(epitope_cleavage_scores, key=lambda x: x[1])
                best_cleavage_position = max_cleavage_score[0]
                best_cleavage_score = max_cleavage_score[1]
                sorted_cleavage_scores = collections.OrderedDict(sorted(epitope_cleavage_scores))
                cleavage_sites = ','.join(['%s:%s' % (key, value) for (key, value) in sorted_cleavage_scores.items()])
            line = current_buffer[sequence_name]
            line.update({
                'Best Cleavage Position': best_cleavage_position,
                'Best Cleavage Score'   : best_cleavage_score,
                'Cleavage Sites'        : cleavage_sites,
            })
            self.write_line(line, writer, current_chunk_id, out_chunk_id, pmp)
        out_chunk_id[0] += 1

    def execute(self):
        chosen_method = str(methods.index(self.method))

        mt_records_dict = pymp.shared.dict(self.get_mt_peptides())

        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter='\t')
            cleavage_cols = ['Best Cleavage Position', 'Best Cleavage Score', 'Cleavage Sites']
            writer = csv.DictWriter(
                output_fh,
                reader.fieldnames if all(col in reader.fieldnames for col in cleavage_cols) else reader.fieldnames+cleavage_cols,
                delimiter='\t',
                lineterminator='\n'
            )
            writer.writeheader()

            # Shared ints
            seq_id = pymp.shared.list([0])
            cycle_ind = pymp.shared.list([1])
            out_chunk_id = pymp.shared.list([0])

            print("Waiting for results from NetChop... |", end='')
            sys.stdout.flush()
            with pymp.Parallel(self.n_threads) as p:
                for chunk_id, chunk in p.iterate(enumerate(lib.run_utils.split_file(reader, 100))): #for chunk in lib.run_utils.split_file(reader, 100):
                    staging_file = tempfile.NamedTemporaryFile(mode='w+')
                    current_buffer = {}
                    seqs_start_diff = {}

                    self.populate_staging_file(chunk, staging_file, current_buffer, seqs_start_diff, mt_records_dict, seq_id)
                    results = self.call_net_chop(staging_file, chosen_method, cycle_ind, p)

                    self.write_results(results, current_buffer, seqs_start_diff, writer, chunk_id, out_chunk_id, p)
                    
            sys.stdout.write('\b\b')
            print("OK")

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            "%s net_chop" % tool,
            description="Predict cleavage sites for neoepitopes.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="Input filtered file with predicted epitopes."
        )
        parser.add_argument(
            'input_fasta',
            help="The required fasta file."
        )
        parser.add_argument(
            'output_file',
            help="Output tsv filename for putative neoepitopes."
        )
        parser.add_argument(
            '--method',
            choices=methods,
            help="NetChop prediction method to use (\"cterm\" for C term 3.0, \"20s\" for 20S 3.0).",
            default='cterm'
        )
        parser.add_argument(
            '--threshold',
            type=float,
            help="NetChop prediction threshold.",
            default=0.5
        )
        return parser


# if __name__ == '__main__':
#     main()
