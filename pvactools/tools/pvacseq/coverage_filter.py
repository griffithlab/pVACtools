import argparse
import sys
import re
import os
import csv

from pvactools.lib.filter import Filter
from pvactools.lib.run_utils import *

def define_parser():
    parser = argparse.ArgumentParser(
        'pvacseq coverage_filter',
        description="Filter variants processed by IEDB by coverage, vaf, and gene expression",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
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
        help="Normal Coverage Cutoff. Sites above this cutoff will be considered.",
        default=5
    )
    parser.add_argument(
        '--tdna-cov', type=int,
        help="Tumor DNA Coverage Cutoff. Sites above this cutoff will be considered.",
        default=10
    )
    parser.add_argument(
        '--trna-cov', type=int,
        help="Tumor RNA Coverage Cutoff. Sites above this cutoff will be considered.",
        default=10
    )
    parser.add_argument(
        '--normal-vaf', type=float_range(0.0,1.0),
        help="Normal VAF Cutoff in decimal format. Sites BELOW this cutoff in normal will be considered.",
        default=0.02
    )
    parser.add_argument(
        '--tdna-vaf', type=float_range(0.0,1.0),
        help="Tumor DNA VAF Cutoff in decimal format. Sites above this cutoff will be considered.",
        default=0.25
    )
    parser.add_argument(
        '--trna-vaf', type=float_range(0.0,1.0),
        help="Tumor RNA VAF Cutoff in decimal format. Sites above this cutoff will be considered.",
        default=0.25
    )
    parser.add_argument(
        '--expn-val', type=float,
        help="Gene and Transcript Expression cutoff. Sites above this cutoff will be considered.",
        default=1.0
    )
    parser.add_argument(
        '--exclude-NAs',
        help="Exclude NA values from the filtered output.",
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
    filter_criteria.append({'column': "Normal Depth", 'operator': '>=', 'threshold': args.normal_cov, 'exclude_nas': args.exclude_NAs})
    filter_criteria.append({'column': "Normal VAF", 'operator': '<=', 'threshold': args.normal_vaf, 'exclude_nas': args.exclude_NAs})
    filter_criteria.append({'column': "Tumor DNA Depth", 'operator': '>=', 'threshold': args.tdna_cov, 'exclude_nas': args.exclude_NAs})
    filter_criteria.append({'column': "Tumor DNA VAF", 'operator': '>=', 'threshold': args.tdna_vaf, 'exclude_nas': args.exclude_NAs})
    filter_criteria.append({'column': "Tumor RNA Depth", 'operator': '>=', 'threshold': args.trna_cov, 'exclude_nas': args.exclude_NAs})
    filter_criteria.append({'column': "Tumor RNA VAF", 'operator': '>=', 'threshold': args.trna_vaf, 'exclude_nas': args.exclude_NAs})
    filter_criteria.append({'column': "Gene Expression", 'operator': '>=', 'threshold': args.expn_val, 'exclude_nas': args.exclude_NAs})
    filter_criteria.append({'column': "Transcript Expression", 'operator': '>=', 'threshold': args.expn_val, 'exclude_nas': args.exclude_NAs})

    Filter(args.input_file, args.output_file, filter_criteria).execute()

if __name__ == "__main__":
    main()
