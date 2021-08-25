import argparse
import sys
import requests
import csv
import tempfile
import re
import os
from time import sleep
import collections
import lib.run_utils
import logging
import random

methods = ['cterm', '20s']

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser("pvacseq net_chop", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'input_file',
        type=argparse.FileType('r'),
        help="Input filtered file with predicted epitopes."
    )
    parser.add_argument(
        'output_file',
        type=argparse.FileType('w'),
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
    args = parser.parse_args(args_input)
    chosen_method = str(methods.index(args.method))
    jobid_searcher = re.compile(r'<!-- jobid: [0-9a-fA-F]*? status: (queued|active)')
    result_delimiter = re.compile(r'-{20,}')
    fail_searcher = re.compile(r'(Failed run|Problematic input:|Unrecognized parameter:)')
    rejected_searcher = re.compile(r'status: rejected')
    success_searcher = re.compile(r'NetChop 3.0 predictions')
    reader = csv.DictReader(args.input_file, delimiter='\t')
    writer = csv.DictWriter(
        args.output_file,
        reader.fieldnames+['Best Cleavage Position', 'Best Cleavage Score', 'Cleavage Sites'],
        delimiter='\t',
        lineterminator='\n'
    )
    writer.writeheader()
    x = 0
    i=1
    for chunk in lib.run_utils.split_file(reader, 100):
        staging_file = tempfile.NamedTemporaryFile(mode='w+')
        current_buffer = {}
        for line in chunk:
            sequence_id = ('%010x'%x)[-10:]
            staging_file.write('>'+sequence_id+'\n')
            if 'Epitope Seq' in line:
                staging_file.write(line['Epitope Seq']+'\n')
            else:
                staging_file.write(line['MT Epitope Seq']+'\n')
            current_buffer[sequence_id] = {k:line[k] for k in line}
            x+=1
        staging_file.seek(0)
        response = query_netchop_server(staging_file, chosen_method, args.threshold, jobid_searcher)

        if fail_searcher.search(response.content.decode()):
            raise Exception("NetChop encountered an error during processing.\n{}".format(response.content.decode()))

        while rejected_searcher.search(response.content.decode()):
            logging.warning("Too many jobs submitted to NetChop server. Waiting to retry.")
            sleep(random.randint(5, 10))
            staging_file.seek(0)
            response = query_netchop_server(staging_file, chosen_method, args.threshold, jobid_searcher)

        if success_searcher.search(response.content.decode()):
            results = [item.strip() for item in result_delimiter.split(response.content.decode())]
            for i in range(2, len(results), 4): #examine only the parts we want, skipping all else
                score = -1
                pos = 0
                sequence_name = False
                cleavage_scores = {}
                for line in results[i].split('\n'):
                    data = [word for word in line.strip().split(' ') if len(word)]
                    if not sequence_name:
                        sequence_name = data[4]
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
                    max_cleavage_score = max(cleavage_scores.items(), key=lambda x: x[1])
                    best_cleavage_position = max_cleavage_score[0]
                    best_cleavage_score = max_cleavage_score[1]
                    sorted_cleavage_scores = collections.OrderedDict(sorted(cleavage_scores.items()))
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
    args.output_file.close()
    args.input_file.close()

def query_netchop_server(staging_file, chosen_method, threshold, jobid_searcher):
    response = requests.post(
        "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi",
        files={'SEQSUB':(staging_file.name, staging_file, 'text/plain')},
        data = {
            'configfile':'/var/www/html/services/NetChop-3.1/webface.cf',
            'SEQPASTE':'',
            'method':chosen_method,
            'thresh':'%0f'%threshold
        }
    )
    if response.status_code != 200:
        raise Exception("Error posting request to NetChop server.\n{}".format(response.content.decode()))
    while jobid_searcher.search(response.content.decode()):
        sleep(10)
        response = requests.get(response.url)
        if response.status_code != 200:
            raise Exception("Error posting request to NetChop server.\n{}".format(response.content.decode()))
    return response

if __name__ == '__main__':
    main()
