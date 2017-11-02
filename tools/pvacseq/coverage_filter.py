import argparse
import sys
import re
import os
import csv
from lib.filter import *

def define_parser():
    parser = argparse.ArgumentParser('pvacseq coverage_filter')
    parser.add_argument(
        'input_file',
        help="The final report .tsv file to filter"
    )
    parser.add_argument(
        'output_file',
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
    parser.add_argument(
        '--exclude-NAs',
        help="Exclude NA values from the filtered output. Default: False",
        default=False,
        action='store_true'
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
    filter_criteria = []
    filter_criteria.append({'column': "Normal_Depth", 'operator': '>=', 'threshold': args.normal_cov})
    filter_criteria.append({'column': "Normal_VAF", 'operator': '<=', 'threshold': args.normal_vaf})
    filter_criteria.append({'column': "Tumor_DNA_Depth", 'operator': '>=', 'threshold': args.tdna_cov})
    filter_criteria.append({'column': "Tumor_DNA_VAF", 'operator': '>=', 'threshold': args.tdna_vaf})
    filter_criteria.append({'column': "Tumor_RNA_Depth", 'operator': '>=', 'threshold': args.trna_cov})
    filter_criteria.append({'column': "Tumor_RNA_VAF", 'operator': '>=', 'threshold': args.trna_vaf})
    filter_criteria.append({'column': "Gene_Expression", 'operator': '>=', 'threshold': args.expn_val})
    filter_criteria.append({'column': "Transcript_Expression", 'operator': '>=', 'threshold': args.expn_val})

    Filter(args.input_file, args.output_file, filter_criteria, args.exclude_NAs).execute()

if __name__ == "__main__":
    main()
