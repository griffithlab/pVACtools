import argparse
import sys
from subprocess import call
import os
import pkg_resources
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
