import argparse
import sys
import re
import os
import csv

from pvactools.lib.filter import Filter, FilterCriterion
from pvactools.lib.run_utils import *

def define_parser():
    parser = argparse.ArgumentParser(
        'pvacsplice coverage_filter',
        description="Filter variants processed by IEDB by coverage, vaf, and gene expression",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'input_file',
        help="The all_epitopes.tsv or filtered.tsv pVACsplice report file to filter."
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
    filter_criteria.append(FilterCriterion("Normal Depth", '>=', args.normal_cov))
    filter_criteria.append(FilterCriterion("Normal VAF", '<=', args.normal_vaf))
    filter_criteria.append(FilterCriterion("Tumor DNA Depth", '>=', args.tdna_cov))
    filter_criteria.append(FilterCriterion("Tumor DNA VAF", '>=', args.tdna_vaf))
    filter_criteria.append(FilterCriterion("Tumor RNA Depth", '>=', args.trna_cov))
    filter_criteria.append(FilterCriterion("Tumor RNA VAF", '>=', args.trna_vaf))
    filter_criteria.append(FilterCriterion("Gene Expression", '>=', args.expn_val))
    filter_criteria.append(FilterCriterion("Transcript Expression", '>=', args.expn_val))

    Filter(args.input_file, args.output_file, filter_criteria).execute()

if __name__ == "__main__":
    main()
