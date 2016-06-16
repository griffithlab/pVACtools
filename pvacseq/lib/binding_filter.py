import argparse
import sys
import re
import os
import csv


def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq binding_filter')
    parser.add_argument('fof', type=argparse.FileType('r'),
                        help="FOF containing list of parsed epitope files " +
                        "for different allele-length combinations (same sample)")
    parser.add_argument('output', type=argparse.FileType('w'),
                        help="Output .xls file containing list of filtered " +
                        "epitopes based on binding affinity for each " +
                        "allele-length combination per gene")
    parser.add_argument('-c', '--fold-change', type=int,
                        help="Minimum fold change between mutant binding " +
                        "score and wild-type score. The default is 0, which " +
                        "filters no results, but 1 is often a sensible " +
                        "default (requiring that binding is better to the MT " +
                        "than WT)",
                        default=0,
                        dest="minimum_fold_change")
    parser.add_argument('-b', '--binding-threshold', type=int,
                        help="Report only epitopes where the mutant allele " +
                        "has ic50 binding scores below this value; default 500",
                        default=500,
                        dest="binding_threshold")

    args = parser.parse_args(args_input)

    prediction = {}
    fieldnames = []

    for line in args.fof:
        filepath = line.rstrip()
        sample = filepath.split(".")[0].replace("_netmhc", "")
        if not os.path.isfile(filepath):
            #if path does not exist, try making it relative to the fof
            filepath = os.path.join(
                os.path.dirname(os.path.abspath(args.fof.name)),
                filepath
            )
        file_handle = open(filepath, mode='r')
        reader = csv.DictReader(file_handle, delimiter='\t')
        if len(fieldnames) == 0:
            fieldnames = reader.fieldnames

        prediction[sample] = {}

        for entry in reader:
            name = entry['Gene Name']
            score = int(entry['MT score'])
            fold_change = 0 if entry['Fold Change'] == 'NA' else float(entry['Fold Change'])

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
        file_handle.close()
    args.fof.close()

    writer = csv.DictWriter(
        args.output,
        fieldnames,
        delimiter = '\t',
        lineterminator = '\n'
    )

    writer.writeheader()

    for entry in (
        #flatten the dictionary structure and iterate over each entry
        item for sample
        in sorted(prediction) for gene
        in sorted(prediction[sample]) for item
        in prediction[sample][gene]['GENES']
    ):
        writer.writerow(entry)

    args.output.close()



if __name__ == "__main__":
    main()
