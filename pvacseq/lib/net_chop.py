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
        'input_file',
        type=argparse.FileType('r'),
        help="Input filtered file with predicted epitopes"
    )
    parser.add_argument(
        'output_file',
        type=argparse.FileType('w'),
        help="Output tsv filename for putative neoepitopes"
    )
    parser.add_argument(
        '--method',
        choices=methods,
        help="NetChop prediction method to use (\"cterm\" for C term 3.0, \"20s\" for 20S 3.0).  Default: \"cterm\" (C term 3.0)",
        default='cterm'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        help="NetChop prediction threshold.  Default: 0.5",
        default=0.5
    )
    args = parser.parse_args(args_input)
    chosen_method = str(methods.index(args.method))
    jobid_searcher = re.compile(r'<!-- jobid: [0-9a-fA-F]*? status: (queued|active)')
    result_delimiter = re.compile(r'-+')
    fail_searcher = re.compile(r'(Failed run|Problematic input:)')
    reader = csv.DictReader(args.input_file, delimiter='\t')
    writer = csv.DictWriter(
        args.output_file,
        reader.fieldnames+['Best Cleavage Position', 'Best Cleavage Score'],
        delimiter='\t',
        lineterminator='\n'
    )
    writer.writeheader()
    x = 0
    i=1
    print("Waiting for results from NetChop... |", end='')
    sys.stdout.flush()
    for chunk in split_file(reader, 100):
        staging_file = tempfile.NamedTemporaryFile(mode='w+')
        current_buffer = {}
        for line in chunk:
            sequence_id = '%010x'%x
            staging_file.write('>'+sequence_id+'\n')
            staging_file.write(line['MT Epitope Seq']+'\n')
            current_buffer[sequence_id] = {k:line[k] for k in line}
            x+=1
        staging_file.seek(0)
        response = requests.post(
            "http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi",
            files={'SEQSUB':(staging_file.name, staging_file, 'text/plain')},
            data = {
                'configfile':'/usr/opt/www/pub/CBS/services/NetChop-3.1/NetChop.cf',
                'SEQPASTE':'',
                'method':chosen_method,
                'thresh':'%0f'%args.threshold
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
            print("NetChop encountered an error during processing")
            sys.exit(1)
        score = -1
        pos = 0
        sequence_name = False
        for line in response.content.decode().split('\n'):
            if result_delimiter.match(line):
                mode = (mode+1)%5
            elif mode==2: #Reading results from the current sequence
                data = [word for word in line.strip().split(' ') if len(word)]
                currentPosition = data[0]
                currentScore = float(data[3])
                if not sequence_name:
                    sequence_name = data[4]
                if currentScore > score:
                    score = currentScore
                    pos = currentPosition
            elif mode==3: #End of results for this sequence.  Output the best sequence
                line = current_buffer[sequence_name]
                line.update({
                    'Best Cleavage Position':pos,
                    'Best Cleavage Score':score
                })
                writer.writerow(line)
                score=-1
                sequence_name=False
                mode+=1
    sys.stdout.write('\b\b')
    print("OK")
    args.output_file.close()
    args.input_file.close()


if __name__ == '__main__':
    main()
