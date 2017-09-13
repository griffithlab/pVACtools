import argparse
import sys
import requests
import csv
import tempfile
import re
import os
from time import sleep

cycle = ['|', '/', '-', '\\']
methods = ['cterm', '20s']

def split_file(reader, lines=400):
    from itertools import islice, chain
    tmp = next(reader)
    while tmp!="":
        yield chain([tmp], islice(reader, lines-1))
        try:
            tmp = next(reader)
        except StopIteration:
            return

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser("pvacseq net_chop")
    parser.add_argument(
        'input_file',
        type=argparse.FileType('r'),
        help="Input filtered file with predicted epitopes"
    )
    parser.add_argument(
        'output_file',
        type=argparse.FileType('w'),
        help="Output TSV filename for putative neoepitopes"
    )
    args = parser.parse_args(args_input)
    jobid_searcher = re.compile(r'<!-- jobid: [0-9a-fA-F]*? status: (queued|active)')
    result_delimiter = re.compile(r'-{20,}')
    fail_searcher = re.compile(r'(Failed run|Problematic input:)')
    allele_searcher = re.compile(r'^(.*?) : Distance to trai?ning data .*? nearest neighbor (.*?)\)$', re.MULTILINE)
    reader = csv.DictReader(args.input_file, delimiter='\t')
    writer = csv.DictWriter(
        args.output_file,
        reader.fieldnames+['Predicted Stability', 'Half Life', 'Stability Rank', 'NetMHCstab allele'],
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
            sequence_id = ('%010x'%x)[-10:]
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
        if fail_searcher.search(response.content.decode()):
            sys.stdout.write('\b\b')
            print('Failed!')
            print("NetMHCStabPan encountered an error during processing")
            sys.exit(1)
        pending = []
        allele_map = {item[0]:item[1] for item in allele_searcher.findall(response.content.decode())}
        results = [item.strip() for item in result_delimiter.split(response.content.decode())]
        for i in range(2, len(results), 4): #examine only the parts we want, skipping all else
            for line in results[i].split('\n'):
                data = [word for word in line.strip().split(' ') if len(word)]
                line = current_buffer[data[3]]
                if data[1] == line['HLA Allele'] and len(data[2]) == int(line['Peptide Length']):
                    line.update({
                        'Predicted Stability':data[4],
                        'Half Life':data[5],
                        'Stability Rank':data[6],
                        'NetMHCstab allele':allele_map[line['HLA Allele'].replace('*', '', 1)]
                    })
                    pending.append([int(data[3], 16), {k:line[k] for k in line}])
        writer.writerows([{k:entry[1][k] for k in entry[1]} for entry in sorted(pending, key=lambda x:x[0])])
    sys.stdout.write('\b\b')
    print("OK")
    args.output_file.close()
    args.input_file.close()


if __name__ == '__main__':
    main()
