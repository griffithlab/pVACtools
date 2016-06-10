import argparse
import csv
import re

parser = argparse.ArgumentParser(description='')
parser.add_argument('input_file', type=argparse.FileType('r'), help="input FASTA file with variant sequences for wildtype(WT) and mutant(MT) proteins generated using 'GenerateVariantSequences.pl' and filtered using 'FilterSeq.pl'")
parser.add_argument('output_file', type=argparse.FileType('w'), help='output Key file for lookup')

args = parser.parse_args()

tsvout = csv.writer(args.output_file, delimiter='\t', lineterminator='\n')

i = 1
pattern = re.compile('>');
for line in args.input_file:
    match = pattern.match(line)
    if match is not None:
        original_name = line.rstrip()
        new_name      = "Entry_%s" % i
        tsvout.writerow([new_name, original_name])
        i += 1

args.input_file.close()
args.output_file.close()
