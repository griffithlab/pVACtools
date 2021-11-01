import argparse
import sys
import requests
import csv
import tempfile
import re
import os
from time import sleep
import lib.run_utils
import random
import logging

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser("pvacseq net_chop", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    fail_searcher = re.compile(r'(Failed run|Problematic input:|Configuration error)')
    rejected_searcher = re.compile(r'status: rejected')
    success_searcher = re.compile(r'Rank Threshold for Strong binding peptides')
    cannot_open_file_searcher = re.compile(r'Cannot open file')
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
    for chunk in lib.run_utils.split_file(reader, 100):
        peptide_lengths = set()
        staging_file = tempfile.NamedTemporaryFile(mode='w+')
        current_buffer = {}
        alleles_in_chunk = set()
        for line in chunk:
            sequence_id = ('%010x'%x)[-10:]
            staging_file.write('>'+sequence_id+'\n')
            if 'Epitope Seq' in line:
                staging_file.write(line['Epitope Seq']+'\n')
                peptide_lengths.add(str(len(line['Epitope Seq'])))
            else:
                staging_file.write(line['MT Epitope Seq']+'\n')
                peptide_lengths.add(str(len(line['MT Epitope Seq'])))
            alleles_in_chunk.add(line['HLA Allele'])
            current_buffer[sequence_id] = {k:line[k] for k in line}
            x+=1
        staging_file.seek(0)
        allele_list = [allele.replace('*', '') for allele in alleles_in_chunk]
        allele_list.sort()
        response = query_netmhcstabpan_server(staging_file, peptide_lengths, allele_list, jobid_searcher)

        if fail_searcher.search(response.content.decode()):
            raise Exception("NetMHCstabpan encountered an error during processing.\n{}".format(response.content.decode()))

        while rejected_searcher.search(response.content.decode()):
            logging.warning("Too many jobs submitted to NetMHCstabpan server. Waiting to retry.")
            sleep(random.randint(5, 10))
            staging_file.seek(0)
            response = query_netmhcstabpan_server(staging_file, peptide_lengths, allele_list, jobid_searcher)

        if cannot_open_file_searcher.search(response.content.decode()):
            sleep(random.randint(5, 10))
            staging_file.seek(0)
            response = query_netmhcstabpan_server(staging_file, peptide_lengths, allele_list, jobid_searcher)
            while rejected_searcher.search(response.content.decode()):
                sleep(random.randint(5, 10))
                staging_file.seek(0)
                response = query_netmhcstabpan_server(staging_file, peptide_lengths, allele_list, jobid_searcher)
            if cannot_open_file_searcher.search(response.content.decode()):
                raise Exception("NetMHCstabpan server was unable to read the submitted fasta file:\n{}.".format(staging_file.read()))

        if success_searcher.search(response.content.decode()):
            pending = []
            allele_map = {item[0]:item[1] for item in allele_searcher.findall(response.content.decode())}
            results = [item.strip() for item in result_delimiter.split(response.content.decode())]
            for i in range(2, len(results), 4): #examine only the parts we want, skipping all else
                for line in results[i].split('\n'):
                    data = [word for word in line.strip().split(' ') if len(word)]
                    line = current_buffer[data[3]]
                    if 'Epitope Seq' in line:
                        length = len(line['Epitope Seq'])
                    else:
                        length = len(line['MT Epitope Seq'])
                    if data[1] == line['HLA Allele'] and len(data[2]) == length:
                        line.update({
                            'Predicted Stability':data[4],
                            'Half Life':data[5],
                            'Stability Rank':data[6],
                            'NetMHCstab allele':allele_map[line['HLA Allele'].replace('*', '', 1)]
                        })
                        pending.append([int(data[3], 16), {k:line[k] for k in line}])
            if len(pending) == 0:
                raise Exception("Unexpected return value from NetMHCstabpan server. Unable to parse response.\n{}".format(response.content.decode()))
            writer.writerows([{k:entry[1][k] for k in entry[1]} for entry in sorted(pending, key=lambda x:x[0])])
        else:
            raise Exception("Unexpected return value from NetMHCstabpan server. Unable to parse response.\n{}".format(response.content.decode()))
    args.output_file.close()
    args.input_file.close()

def query_netmhcstabpan_server(staging_file, peptide_lengths, allele_list, jobid_searcher):
    response = requests.post(
        "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi",
        files={'SEQSUB':(staging_file.name, staging_file, 'text/plain')},
        data = {
            'configfile':'/var/www/html/services/NetMHCstabpan-1.0/webface.cf',
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
            'sort2':'-1',
        }
    )
    if response.status_code != 200:
        raise Exception("Error posting request to NetMHCstabpan server.\n{}".format(response.content.decode()))
    while jobid_searcher.search(response.content.decode()):
        sleep(10)
        response = requests.get(response.url)
        if response.status_code != 200:
            raise Exception("Error posting request to NetMHCstabpan server.\n{}".format(response.content.decode()))
    return response


if __name__ == '__main__':
    main()
