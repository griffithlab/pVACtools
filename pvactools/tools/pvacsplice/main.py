import argparse
import sys
from subprocess import call
import os
from pvactools.tools.pvacsplice import *

def define_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    #add subcommands
    run_main_program_parser = subparsers.add_parser(
        "run",
        help="Run the pVACsplice pipeline.",
        add_help=False
    )
    run_main_program_parser.set_defaults(func=run)

    binding_filter_parser = subparsers.add_parser(
        "binding_filter",
        help="Filter variants processed by IEDB by binding score.",
        add_help=False
    )
    binding_filter_parser.set_defaults(func=binding_filter)

    coverage_filter_parser = subparsers.add_parser(
        "coverage_filter",
        help="Filter variants processed by IEDB by coverage, vaf, and gene expression.",
        add_help=False
    )
    coverage_filter_parser.set_defaults(func=coverage_filter)

    transcript_support_level_filter_parser = subparsers.add_parser(
        "transcript_support_level_filter",
        help="Filter variants processed by IEDB by transcript support level.",
        add_help=False
    )
    transcript_support_level_filter_parser.set_defaults(func=transcript_support_level_filter)

    top_score_filter_parser = subparsers.add_parser(
        "top_score_filter",
        help="Pick the best neoepitope for each variant.",
        add_help=False,
    )
    top_score_filter_parser.set_defaults(func=top_score_filter)

    gtf_to_tsv = subparsers.add_parser(
        "gtf_to_tsv",
        help="Filter and convert reference GTF to TSV format.",
        add_help=False
    )
    gtf_to_tsv.set_defaults(func=run)
    
    return parser


def main():
    parser = define_parser()
    args = parser.parse_known_args()
    try:
        args[0].func.main(args[1])
    except AttributeError as e:
        parser.print_help()
        print("Error: No command specified")
        sys.exit(-1)


if __name__ == '__main__':
    main()
