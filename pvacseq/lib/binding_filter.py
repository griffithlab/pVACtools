import argparse
import sys
import re
import os
import csv


def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq binding_filter')
    parser.add_argument('input_files', type=argparse.FileType('r'),
                        nargs="+",
                        help="List of parsed epitope files " +
                        "for different allele-length combinations (same sample)")
    parser.add_argument('output_file', type=argparse.FileType('w'),
                        help="Output .xls file containing list of filtered " +
                        "epitopes based on binding affinity for each " +
                        "allele-length combination per gene")
    parser.add_argument('-c', '--minimum-fold-change', type=int,
                        help="Minimum fold change between mutant binding " +
                        "score and wild-type score. The default is 0, which " +
                        "filters no results, but 1 is often a sensible " +
                        "default (requiring that binding is better to the MT " +
                        "than WT)",
                        default=0)
    parser.add_argument('-b', '--binding-threshold', type=int,
                        help="Report only epitopes where the mutant allele " +
                        "has ic50 binding scores below this value; default 500",
                        default=500)
    parser.add_argument('-m', '--top-score-metric',
                        choices=['lowest', 'median'],
                        default='median',
                        help="Which ic50 scoring metric to use when filtering epitopes by binding-threshold. " +
                        "lowest: Best MT Score - lowest WT ic50 binding score between all chosen prediction methods. " +
                        "median: Median MT Score - median WT ic50 binding score between all chosen prediction methods. " +
                        "Default: median")

    args = parser.parse_args(args_input)

    prediction = {}
    fieldnames = []

    for input_file in args.input_files:
        sample = os.path.basename(input_file.name).split(".")[0].replace("_netmhc", "")
        reader = csv.DictReader(input_file, delimiter='\t')
        if len(fieldnames) == 0:
            fieldnames = reader.fieldnames

        if sample not in prediction:
            prediction[sample] = {}

        for entry in reader:
            name = entry['Gene Name']
            if args.top_score_metric == 'median':
                score = float(entry['Median MT Score All Methods'])
            elif args.top_score_metric == 'best':
                score = float(entry['Best MT Score'])
            fold_change = sys.maxsize if entry['Fold Change'] == 'NA' else float(entry['Fold Change'])

            if score > args.binding_threshold or fold_change < args.minimum_fold_change:
                continue

            if (name not in prediction[sample] or
                    score < prediction[sample][name]['SCORE']):
                prediction[sample][name] = {
                    'GENES' : [entry],
                    'SCORE' : score
                }
            elif score == prediction[sample][name]['SCORE']:
                prediction[sample][name]['GENES'].append(entry)
        input_file.close()

    writer = csv.DictWriter(
        args.output_file,
        fieldnames,
        delimiter = '\t',
        lineterminator = '\n'
    )

    writer.writeheader()

    writer.writerows(
        entry
        for sample in sorted(prediction)
        for gene in sorted(prediction[sample])
        for entry in prediction[sample][gene]['GENES']
    )

    args.output_file.close()



if __name__ == "__main__":
    main()
