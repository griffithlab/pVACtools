import argparse
import sys
import re
import os
import csv

def define_parser():
    parser = argparse.ArgumentParser('pvacseq coverage_filter')
    parser.add_argument(
        'input_file', type=argparse.FileType('r'),
        help="The final report .tsv file to filter"
    )
    parser.add_argument(
        'output_file', type=argparse.FileType('w'),
        help="Output .tsv file containing list of filtered epitopes based on coverage and expression values"
    )
    parser.add_argument(
        '--normal-cov', type=int,
        help="Normal Coverage Cutoff. Sites above this cutoff will be considered. "
             + "Default: 5",
        default=5
    )
    parser.add_argument(
        '--tdna-cov', type=int,
        help="Tumor DNA Coverage Cutoff. Sites above this cutoff will be considered. "
             +"Default: 10",
        default=10
    )
    parser.add_argument(
        '--trna-cov', type=int,
        help="Tumor RNA Coverage Cutoff. Sites above this cutoff will be considered. "
             + "Default: 10",
        default=10
    )
    parser.add_argument(
        '--normal-vaf', type=int,
        help="Normal VAF Cutoff. Sites BELOW this cutoff in normal will be considered. "
             + "Default: 2",
        default=2
    )
    parser.add_argument(
        '--tdna-vaf', type=int,
        help="Tumor DNA VAF Cutoff. Sites above this cutoff will be considered. "
             + "Default: 40",
        default=40
    )
    parser.add_argument(
        '--trna-vaf', type=int,
        help="Tumor RNA VAF Cutoff. Sites above this cutoff will be considered. "
             + "Default: 40",
        default=40
    )
    parser.add_argument(
        '--expn-val', type=int,
        help="Gene and Transcript Expression cutoff. Sites above this cutoff will be considered"
             + "Default: 1",
         default=1
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

#### COVERAGE COLUMNS ##
#Normal Depth
#Normal VAF
#Tumor DNA Depth
#Tumor DNA VAF
#Tumor RNA Depth
#Tumor RNA VAF
#Gene Expression
#Transcript Expression

    reader = csv.DictReader(args.input_file, delimiter='\t')
    writer = csv.DictWriter(
        args.output_file,
        reader.fieldnames,
        delimiter = '\t',
        lineterminator = '\n'
    )
    writer.writeheader()
    for entry in reader:
        if ('Normal Depth' in entry
            and entry['Normal Depth']
            and entry['Normal Depth'] != 'NA'
            and float(entry ['Normal Depth']) < args.normal_cov):
            continue

        if ('Normal VAF' in entry
            and entry['Normal VAF']
            and entry['Normal VAF'] != 'NA'
            and float(entry['Normal VAF']) > args.normal_vaf):
            continue

        if ('Tumor DNA Depth' in entry
            and entry['Tumor DNA Depth']
            and entry['Tumor DNA Depth'] != 'NA'
            and float(entry['Tumor DNA Depth']) < args.tdna_cov):
            continue

        if ('Tumor DNA VAF' in entry
            and entry['Tumor DNA VAF']
            and entry['Tumor DNA VAF'] != 'NA'
            and float(entry['Tumor DNA VAF']) < args.tdna_vaf):
            continue

        if ('Tumor RNA Depth' in entry
            and entry['Tumor RNA Depth']
            and entry['Tumor RNA Depth'] != 'NA'
            and float(entry['Tumor RNA Depth']) < args.trna_cov):
            continue

        if ('Tumor RNA VAF' in entry
            and entry['Tumor RNA VAF']
            and entry['Tumor RNA VAF'] != 'NA'
            and float(entry['Tumor RNA VAF']) < args.trna_vaf):
            continue

        if ('Gene Expression' in entry
            and entry['Gene Expression']
            and entry['Gene Expression'] != 'NA'
            and float(entry['Gene Expression']) < args.expn_val):
            continue

        if ('Transcript Expression' in entry
            and entry['Transcript Expression']
            and entry['Transcript Expression'] != 'NA'
            and float(entry['Transcript Expression']) < args.expn_val):
            continue

        writer.writerow(entry)

    args.input_file.close()
    args.output_file.close()

if __name__ == "__main__":
    main()
