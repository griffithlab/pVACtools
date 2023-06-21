import argparse
import sys
import re
import os
import csv

from pvactools.lib.filter import Filter, FilterCriterion
from pvactools.lib.run_utils import *

def define_parser():
    parser = argparse.ArgumentParser(
        'pvacfuse coverage_filter',
        description="Filter variants processed by IEDB by read support and expression",
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
        '--read-support', type=int,
        help="Read Support Cutoff. Sites above this cutoff will be considered.",
        default=5
    )
    parser.add_argument(
        '--expn-val', type=float,
        help="Expression Cutoff. Expression is meassured as FFPM (fusion fragments per million total reads). Sites above this cutoff will be considered.",
        default=0.1
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
#Read Support
#Expression
    filter_criteria = []
    filter_criteria.append(FilterCriterion("Read Support", '>=', args.read_support, exclude_nas=args.exclude_NAs))
    filter_criteria.append(FilterCriterion("Expression", '>=', args.expn_val, exclude_nas=args.exclude_NAs))

    Filter(args.input_file, args.output_file, filter_criteria).execute()

if __name__ == "__main__":
    main()
