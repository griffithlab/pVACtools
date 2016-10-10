import argparse
import sys
import requests
import csv
import tempfile
import re
import os
from time import sleep
import collections

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
    result_delimiter = re.compile(r'-{20,}')
    fail_searcher = re.compile(r'(Failed run|Problematic input:)')
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
    print("Waiting for results from NetChop... |", end='')
    sys.stdout.flush()
    for chunk in split_file(reader, 100):
        staging_file = tempfile.NamedTemporaryFile(mode='w+')
        current_buffer = {}
        for line in chunk:
            sequence_id = ('%010x'%x)[-10:]
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

        results = [item.strip() for item in result_delimiter.split(response.content.decode())]
        for i in range(2, len(results), 4): #examine only the parts we want, skipping all else
            score = -1
            pos = 0
            sequence_name = False
            cleavage_scores = {}
            for line in results[i].split('\n'):
                data = [word for word in line.strip().split(' ') if len(word)]
                currentPosition = data[0]
                isCleavage = data[2]
                if isCleavage is not 'S':
                    continue
                currentScore = float(data[3])
                if not sequence_name:
                    sequence_name = data[4]
                cleavage_scores[currentPosition] = currentScore
            max_cleavage_score = max(cleavage_scores.items(), key=lambda x: x[1])
            sorted_cleavage_scores = collections.OrderedDict(sorted(cleavage_scores.items()))
            line = current_buffer[sequence_name]
            line.update({
                'Best Cleavage Position':max_cleavage_score[0],
                'Best Cleavage Score'   :max_cleavage_score[1],
                'Cleavage Sites'        :','.join(['%s:%s' % (key, value) for (key, value) in sorted_cleavage_scores.items()])
            })
            writer.writerow(line)
    sys.stdout.write('\b\b')
    print("OK")
    args.output_file.close()
    args.input_file.close()


if __name__ == '__main__':
    main()
