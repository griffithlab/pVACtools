import argparse
import sys
import requests
import csv
import tempfile
import re
import os
from time import sleep
from lib.main import split_file

cycle = ['|', '/', '-', '\\']
methods = ['cterm', '20s']

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser("pvacseq net_chop")
    parser.add_argument(
        'input',
        type=argparse.FileType('r'),
        help="Input filtered file with predicted epitopes"
    )
    parser.add_argument(
        'output',
        type=argparse.FileType('w'),
        help="Output FASTA filename for putative neoepitopes"
    )
    args = parser.parse_args(args_input)
    jobid_searcher = re.compile(r'<!-- jobid: [0-9a-fA-F]*? status: (queued|active)')
    result_delimiter = re.compile(r'-+')
    fail_searcher = re.compile(r'(Failed run|Problematic input:)')
    reader = csv.DictReader(args.input, delimiter='\t')
    writer = csv.DictWriter(
        args.output,
        reader.fieldnames+['Predicted Stability', 'Half Life', 'Stability Rank'],
        delimiter='\t',
        lineterminator='\n'
    )
    writer.writeheader()
    x = 0
    i=1
    print("Waiting for results from NetMHCStabPan... |", end='')
    sys.stdout.flush()
    for chunk in split_file(reader, 100):
        peptide_lengths = set()
        staging_file = tempfile.NamedTemporaryFile(mode='w+')
        current_buffer = {}
        alleles_in_chunk = set()
        for line in chunk:
            sequence_id = '%010x'%x
            staging_file.write('>'+sequence_id+'\n')
            staging_file.write(line['MT Epitope Seq']+'\n')
            alleles_in_chunk.add(line['HLA Allele'])
            peptide_lengths.add(line['Peptide Length'])
            current_buffer[sequence_id] = {k:line[k] for k in line}
            x+=1
        staging_file.seek(0)
        allele_list = [allele.replace('*', '') for allele in alleles_in_chunk]
        allele_list.sort()
        response = requests.post(
            "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi",
            files={'SEQSUB':(staging_file.name, staging_file, 'text/plain')},
            data = {
                'configfile':'/usr/opt/www/pub/CBS/services/NetMHCstabpan-1.0/NetMHCstabpan.cf',
                'inp':'0',
                'len': ','.join(peptide_lengths),
                'master':'1',
                'slave0':allele_list[-1],
                'allele':','.join(allele_list),
                'thrs':'0.5',
                'thrw': '2',
                'incaff': '0',
                'sort1':'-1',
                'waff':'0.8',
                'sort2':'-1'
            }
        )
        while jobid_searcher.search(response.content.decode()):
            for _ in range(10):
                sys.stdout.write('\b'+cycle[i%4])
                sys.stdout.flush()
                sleep(1)
                i+=1
            response = requests.get(response.url)
        mode=0
        if fail_searcher.search(response.content.decode()):
            sys.stdout.write('\b\b')
            print('Failed!')
            print("NetMHCStabPan encountered an error during processing")
            sys.exit(1)
        pending = []
        for line in response.content.decode().split('\n'):
            if result_delimiter.match(line):
                mode = (mode+1)%4
            elif mode==2: #Reading results from the current sequence
                data = [word for word in line.strip().split(' ') if len(word)]
                line = current_buffer[data[3]]
                if data[1] == line['HLA Allele'] and len(data[2]) == int(line['Peptide Length']):
                    line.update({
                        'Predicted Stability':data[4],
                        'Half Life':data[5],
                        'Stability Rank':data[6]
                    })
                    pending.append([int(data[3]), {k:line[k] for k in line}])
        writer.writerows([{k:entry[1][k] for k in entry[1]} for entry in sorted(pending, key=lambda x:x[0])])
    sys.stdout.write('\b\b')
    print("OK")
    args.output.close()
    args.input.close()


if __name__ == '__main__':
    main()
