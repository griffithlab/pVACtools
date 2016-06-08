import argparse
import sys
from subprocess import call
import os
try:
    from . import lib
except SystemError:
    import lib

def main():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()

    #add subcommands
    variant_sequences_parser = subparsers.add_parser("generate_fasta",
                                                     help="Generates a variant peptide FASTA file from the TSV input file",
                                                     add_help=False)
    variant_sequences_parser.set_defaults(func=lib.generate_fasta)

    binding_filter_parser = subparsers.add_parser("binding_filter",
                                                  help="Filters variants processed by NetMHC",
                                                  add_help=False)
    binding_filter_parser.set_defaults(func=lib.binding_filter)

    fasta_key_parser = subparsers.add_parser("generate_fasta_key",
                                             help="Generates a FASTA key file",
                                             add_help=False)
    fasta_key_parser.set_defaults(func=lib.generate_fasta_key)

    parse_netmhc_parser = subparsers.add_parser("parse_output",
                                                help="Parses output from NetMHC",
                                                add_help=False)
    parse_netmhc_parser.set_defaults(func=lib.parse_output)

    run_main_program_parser = subparsers.add_parser("run",
                                                    help="Runs the pVAC-Seq pipeline",
                                                    add_help=False)
    run_main_program_parser.set_defaults(func=lib.main)

    args = parser.parse_known_args()
    try:
        args[0].func.main(args[1])
    except AttributeError as e:
        parser.print_help()
        print("Error: No command specified")
        sys.exit(-1)


if __name__ == '__main__':
    main()
