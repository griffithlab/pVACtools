import argparse
import sys
from tools.pvacvector import *

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    #add subcommands
    run_main_program_parser = subparsers.add_parser(
        "run",
        help="Runs the pVACvector pipeline",
        add_help=False
    )
    run_main_program_parser.set_defaults(func=run)

    valid_alleles_parser = subparsers.add_parser(
        "valid_alleles",
        help="Shows a list of valid allele names",
        add_help=False
    )
    valid_alleles_parser.set_defaults(func=valid_alleles)

    download_example_data_parser = subparsers.add_parser(
        "download_example_data",
        help="Downloads example input and output files",
        add_help=False
    )
    download_example_data_parser.set_defaults(func=download_example_data)

    args = parser.parse_known_args()
    try:
        args[0].func.main(args[1])
    except AttributeError as e:
        parser.print_help()
        print("Error: No command specified")
        sys.exit(-1)

if __name__ == '__main__':
    main()
