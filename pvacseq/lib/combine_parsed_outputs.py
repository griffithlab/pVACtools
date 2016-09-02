import argparse
import sys
import csv

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq combine_parsed_outputs')
    parser.add_argument('input_files', type=argparse.FileType('r'),
                        nargs="+",
                        help="List of parsed epitope files " +
                        "for different allele-length combinations (same sample)")
    parser.add_argument('output_file', type=argparse.FileType('w'),
                        help="Combined output .tsv file")
    args = parser.parse_args(args_input)

    fieldnames = []
    for input_file in args.input_files:
        reader = csv.DictReader(input_file, delimiter='\t')
        if len(fieldnames) == 0:
            fieldnames = reader.fieldnames
        else:
            for fieldname in reader.fieldnames:
                if fieldname not in fieldnames:
                    fieldnames.append(fieldname)
        input_file.seek(0)

    tsv_writer = csv.DictWriter(args.output_file, list(fieldnames), delimiter = '\t', lineterminator = '\n')
    tsv_writer.writeheader()
    for input_file in args.input_files:
        reader = csv.DictReader(input_file, delimiter='\t')
        for row in reader:
            for fieldname in fieldnames:
                if fieldname not in row:
                    row[fieldname] = 'NA'
            tsv_writer.writerow(row)
        input_file.close()

    args.output_file.close()

if __name__ == "__main__":
    main()
