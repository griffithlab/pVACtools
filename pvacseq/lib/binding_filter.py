import argparse
import sys
import re
import csv

def define_parser():
    parser = argparse.ArgumentParser('pvacseq binding_filter')
    parser.add_argument(
        'input_file', type=argparse.FileType('r'),
        help="The final report .tsv file to filter"
    )
    parser.add_argument(
        'output_file', type=argparse.FileType('w'),
        help="Output .tsv file containing list of filtered "
             + "epitopes based on binding affinity"
    )
    parser.add_argument(
        '-b', '--binding-threshold', type=int,
        help="Report only epitopes where the mutant allele "
             + "has ic50 binding scores below this value. Default: 500",
        default=500
    )
    parser.add_argument(
        '-c', '--minimum-fold-change', type=int,
        help="Minimum fold change between mutant binding "
             + "score and wild-type score. The default is 0, which "
             + "filters no results, but 1 is often a sensible "
             + "option (requiring that binding is better to the MT than WT). "
             + "Default: 0",
        default=0
    )
    parser.add_argument(
        '-m', '--top-score-metric',
        choices=['lowest', 'median'],
        help="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. "
             + "lowest: Best MT Score/Corresponding Fold Change - lowest MT ic50 binding score/corresponding fold change of all chosen prediction methods. "
             + "median: Median MT Score/Median Fold Change - median MT ic50 binding score/fold change of all chosen prediction methods. "
             + "Default: median",
        default='median',
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    reader = csv.DictReader(args.input_file, delimiter='\t')
    fieldnames = reader.fieldnames

    writer = csv.DictWriter(
        args.output_file,
        fieldnames,
        delimiter = '\t',
        lineterminator = '\n'
    )
    writer.writeheader()

    for entry in reader:
        name = entry['Gene Name']
        if args.top_score_metric == 'median':
            score = float(entry['Median MT Score'])
            fold_change = sys.maxsize if entry['Median Fold Change'] == 'NA' else float(entry['Median Fold Change'])
        elif args.top_score_metric == 'lowest':
            score = float(entry['Best MT Score'])
            fold_change = sys.maxsize if entry['Corresponding Fold Change'] == 'NA' else float(entry['Corresponding Fold Change'])

        if score > args.binding_threshold or fold_change < args.minimum_fold_change:
            continue

        writer.writerow(entry)

    args.input_file.close()
    args.output_file.close()

if __name__ == "__main__":
    main()
