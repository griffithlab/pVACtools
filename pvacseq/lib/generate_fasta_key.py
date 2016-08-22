import argparse
import csv
import re
import sys
import os

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq generate_fasta_key')
    parser.add_argument('input_file', type=argparse.FileType('r'), help="input FASTA file with variant sequences for wildtype(WT) and mutant(MT) proteins generated using 'GenerateVariantSequences.pl' and filtered using 'FilterSeq.pl'")
    parser.add_argument('output_file', help='output Key file for lookup')

    args = parser.parse_args(args_input)

    tmp_output_file = args.output_file + '.tmp'
    tmp_output_filehandle = open(tmp_output_file, 'w')
    tsvout = csv.writer(tmp_output_filehandle, delimiter='\t', lineterminator='\n')

    i = 1
    pattern = re.compile('>')
    for line in args.input_file:
        match = pattern.match(line)
        if match is not None:
            original_name = line.rstrip()
            new_name      = i
            tsvout.writerow([new_name, original_name])
            i += 1

    tmp_output_filehandle.close()
    os.replace(tmp_output_file, args.output_file)

    args.input_file.close()

if __name__ == '__main__':
    main()
