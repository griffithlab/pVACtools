import argparse
import sys
import csv
import pvactools.lib.sort

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq combine_parsed_outputs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'input_files',
        nargs="+",
        help="List of parsed epitope files for different allele-length combinations (same sample)."
    )
    parser.add_argument(
        'output_file', type=argparse.FileType('w'),
        help="Combined output .tsv file."
    )
    parser.add_argument(
        '--top-score-metric',
        choices=['lowest', 'median'],
        default='median',
        help="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. "
             + "lowest: Use the best MT Score and Corresponding Fold Change (i.e. the lowest MT ic50 binding score and corresponding fold change of all chosen prediction methods). "
             + "median: Use the median MT Score and Median Fold Change (i.e. the median MT ic50 binding score and fold change of all chosen prediction methods).",
    )
    parser.add_argument(
        '--file-type',
        choices=['pVACseq', 'pVACfuse', 'pVACbind'],
        default='pVACseq',
        help="Pipeline that created files to be combined."
    )
    args = parser.parse_args(args_input)

    fieldnames = []
    for input_file in args.input_files:
        with open(input_file, 'r') as input_file_handle:
            reader = csv.DictReader(input_file_handle, delimiter='\t')
            if len(fieldnames) == 0:
                fieldnames = reader.fieldnames
            else:
                for fieldname in reader.fieldnames:
                    if fieldname not in fieldnames:
                        fieldnames.append(fieldname)

    tsv_writer = csv.DictWriter(args.output_file, list(fieldnames), delimiter = '\t', lineterminator = '\n', restval='NA')
    tsv_writer.writeheader()
    for input_file in args.input_files:
        with open(input_file, 'r') as input_file_handle:
            reader = csv.DictReader(input_file_handle, delimiter='\t')
            for row in reader:
                tsv_writer.writerow(row)
    args.output_file.close()

if __name__ == "__main__":
    main()
