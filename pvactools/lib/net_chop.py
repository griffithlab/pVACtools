import argparse
import sys
import requests
from requests.exceptions import Timeout
import csv
import tempfile
import re
import os
from time import sleep
from Bio import SeqIO
import collections
import logging
import random
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

import pvactools.lib.run_utils

methods = ['cterm', '20s']

class NetChop:
    def __init__(self, input_file, input_fasta, output_file, method='cterm', threshold=0.5, file_type='pVACseq'):
        self.input_file = input_file
        self.input_fasta = input_fasta
        self.output_file = output_file
        self.method = method
        self.threshold = float(threshold)
        self.file_type = file_type

    def get_mt_peptides(self):
        records = list(SeqIO.parse(self.input_fasta, "fasta"))
        if self.file_type == 'pVACseq':
            records_dict = {x.id.replace('MT.', ''): str(x.seq) for x in filter(lambda x: x.id.startswith('MT.'), records)}
        else:
            records_dict = {x.id: str(x.seq) for x in records}
        return records_dict

    def extract_flanked_epitope(self, full_peptide, epitope, seq_id):
        flanking_sequence_length = 9
        if epitope not in full_peptide:
            raise Exception("FASTA entry {} ({}) does not contain epitope {}. Please check that the FASTA file matches the input TSV.".format(seq_id, full_peptide, epitope))
        ep_start = full_peptide.index(epitope)
        start = ep_start - flanking_sequence_length
        if start < 0:
            start = 0
        start_diff = ep_start - start
        end = ep_start + len(epitope) + flanking_sequence_length
        return full_peptide[start:end], start_diff

    def execute(self):
        chosen_method = str(methods.index(self.method))
        jobid_searcher = re.compile(r'<!-- jobid: [0-9a-fA-F]*? status: (queued|active)')
        result_delimiter = re.compile(r'-{20,}')

        fail_searcher = re.compile(r'(Failed run|Problematic input:|Unrecognized parameter:)')
        rejected_searcher = re.compile(r'status: rejected')
        success_searcher = re.compile(r'NetChop 3.0 predictions')

        mt_records_dict = self.get_mt_peptides()
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
            x = 0
            i=1
            http = self.setup_adapter()
            for chunk in pvactools.lib.run_utils.split_file(reader, 100):
                staging_file = tempfile.NamedTemporaryFile(mode='w+')
                current_buffer = {}
                seqs_start_diff = {}
                for line in chunk:
                    sequence_id = ('%010x'%x)[-10:]
                    staging_file.write('>'+sequence_id+'\n')
                    if self.file_type == 'pVACbind' or self.file_type == 'pVACfuse':
                        index = line['Mutation']
                        epitope = line['Epitope Seq']
                    else:
                        index = line['Index']
                        epitope = line['MT Epitope Seq']
                    if index not in mt_records_dict:
                        raise Exception("FASTA entry for index {} not found. Please check that the FASTA file matches the input TSV.".format(index))
                    full_peptide = mt_records_dict[index]
                    peptide, start_diff = self.extract_flanked_epitope(full_peptide, epitope, index)
                    staging_file.write(peptide+'\n')
                    current_buffer[sequence_id] = {k:line[k] for k in line}
                    seqs_start_diff[sequence_id] = (start_diff, len(epitope))
                    x+=1
                staging_file.seek(0)
                response = self.query_netchop_server(http, staging_file, chosen_method, self.threshold, jobid_searcher)

                if fail_searcher.search(response.content.decode()):
                    raise Exception("NetChop encountered an error during processing.\n{}".format(response.content.decode()))

                while rejected_searcher.search(response.content.decode()):
                    logging.warning("Too many jobs submitted to NetChop server. Waiting to retry.")
                    sleep(random.randint(5, 10))
                    staging_file.seek(0)
                    response = self.query_netchop_server(http, staging_file, chosen_method, self.threshold, jobid_searcher)

                if success_searcher.search(response.content.decode()):
                    results = [item.strip() for item in result_delimiter.split(response.content.decode())]
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
                            if isCleavage != 'S':
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
                            if len(epitope_cleavage_scores) == 0:
                                best_cleavage_position = 'NA'
                                best_cleavage_score = 'NA'
                                cleavage_sites = 'NA'
                            else:
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
                        writer.writerow(line)
                else:
                    raise Exception("Unexpected return value from NetChop server. Unable to parse response.\n{}".format(response.content.decode()))
            http.close()

    def setup_adapter(self):
        retry_strategy = Retry(
            total=3,
            status_forcelist=[408, 429, 500, 502, 503, 504],
            allowed_methods=["POST", "GET"]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        http = requests.Session()
        http.mount("https://", adapter)
        http.mount("http://", adapter)
        return http

    def query_netchop_server(self, http, staging_file, chosen_method, threshold, jobid_searcher):
        try:
            response = self.post_query(http, staging_file, chosen_method, threshold)
        except Timeout:
            raise Exception("Timeout while posting request to NetChop server. The server may be unresponsive. Please try again later.")
        if response.status_code != 200:
            raise Exception("Error posting request to NetChop server.\n{}".format(response.content.decode()))
        while jobid_searcher.search(response.content.decode()):
            sleep(10)
            try:
                response = http.get(response.url, timeout=(10,60))
            except Timeout:
                raise Exception("Timeout while posting request to NetChop server. The server may be unresponsive. Please try again later.")
            if response.status_code != 200:
                raise Exception("Error posting request to NetChop server.\n{}".format(response.content.decode()))
        return response

    def post_query(self, http, staging_file, chosen_method, threshold):
        response = http.post(
            "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi",
            files={'SEQSUB':(staging_file.name, staging_file, 'text/plain')},
            data = {
                'configfile':'/var/www/html/services/NetChop-3.1/webface.cf',
                'SEQPASTE':'',
                'method':chosen_method,
                'thresh':'%0f'%threshold
            },
            timeout=(10,60)
        )
        return response

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
