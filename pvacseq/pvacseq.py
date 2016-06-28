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
    run_main_program_parser = subparsers.add_parser("run",
                                                    help="Runs the pVAC-Seq pipeline",
                                                    add_help=False)
    run_main_program_parser.set_defaults(func=lib.main)

    convert_vcf_parser = subparsers.add_parser("convert_vcf",
                                               help="Converts a VCF into TSV format for downstream steps",
                                               add_help=False)
    convert_vcf_parser.set_defaults(func=lib.convert_vcf)

    variant_sequences_parser = subparsers.add_parser("generate_fasta",
                                                     help="Generates a variant peptide FASTA file from the TSV input file",
                                                     add_help=False)
    variant_sequences_parser.set_defaults(func=lib.generate_fasta)

    fasta_key_parser = subparsers.add_parser("generate_fasta_key",
                                             help="Generates a FASTA key file",
                                             add_help=False)
    fasta_key_parser.set_defaults(func=lib.generate_fasta_key)

    parse_netmhc_parser = subparsers.add_parser("parse_output",
                                                help="Parses output from NetMHC",
                                                add_help=False)
    parse_netmhc_parser.set_defaults(func=lib.parse_output)

    binding_filter_parser = subparsers.add_parser("binding_filter",
                                                  help="Filters variants processed by NetMHC by binding score",
                                                  add_help=False)
    binding_filter_parser.set_defaults(func=lib.binding_filter)

    coverage_filter_parser = subparsers.add_parser("coverage_filter",
                                                   help="Filters variants processed by NetMHC by coverage, vaf, and gene expression",
                                                   add_help=False)
    coverage_filter_parser.set_defaults(func=lib.coverage_filter)

    download_example_data_parser = subparsers.add_parser("download_example_data",
                                                         help="Downloads example input and output files",
                                                         add_help=False)
    download_example_data_parser.set_defaults(func=lib.download_example_data)

    install_vep_plugin_parser = subparsers.add_parser("install_vep_plugin",
                                                      help="Installs the Wildtype VEP plugin into your VEP_plugins directory",
                                                      add_help=False)
    install_vep_plugin_parser.set_defaults(func=lib.install_vep_plugin)

    args = parser.parse_known_args()
    try:
        args[0].func.main(args[1])
    except AttributeError as e:
        parser.print_help()
        print("Error: No command specified")
        sys.exit(-1)


if __name__ == '__main__':
    main()
