import argparse
import sys
from subprocess import call
import os
import pkg_resources
from tools.pvacseq import *

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    #add subcommands
    run_main_program_parser = subparsers.add_parser(
        "run",
        help="Runs the pVACseq pipeline",
        add_help=False
    )
    run_main_program_parser.set_defaults(func=run)

    binding_filter_parser = subparsers.add_parser(
        "binding_filter",
        help="Filters variants processed by IEDB by binding score",
        add_help=False
    )
    binding_filter_parser.set_defaults(func=binding_filter)

    coverage_filter_parser = subparsers.add_parser(
        "coverage_filter",
        help="Filters variants processed by IEDB by coverage, vaf, and gene expression",
        add_help=False
    )
    coverage_filter_parser.set_defaults(func=coverage_filter)

    top_score_filter_parser = subparsers.add_parser(
        "top_score_filter",
        help="Pick the best neoepitope for each variant",
        add_help=False,
    )
    top_score_filter_parser.set_defaults(func=top_score_filter)

    generate_protein_fasta_parser = subparsers.add_parser(
        "generate_protein_fasta",
        help="Generate an annotated fasta file from a VCF with protein sequences of mutations and matching wildtypes",
        add_help=False
    )
    generate_protein_fasta_parser.set_defaults(func=generate_protein_fasta)

    download_example_data_parser = subparsers.add_parser(
        "download_example_data",
        help="Downloads example input and output files",
        add_help=False
    )
    download_example_data_parser.set_defaults(func=download_example_data)

    install_vep_plugin_parser = subparsers.add_parser(
        "install_vep_plugin",
        help="Installs the Wildtype VEP plugin into your VEP_plugins directory",
        add_help=False
    )
    install_vep_plugin_parser.set_defaults(func=install_vep_plugin)

    valid_alleles_parser = subparsers.add_parser(
        "valid_alleles",
        help="Shows a list of valid allele names",
        add_help=False
    )
    valid_alleles_parser.set_defaults(func=valid_alleles)

    config_files_parser = subparsers.add_parser(
        "config_files",
        help="Documentation for the configuration files",
        add_help=False
    )
    config_files_parser.set_defaults(func=config_files)

    args = parser.parse_known_args()
    try:
        args[0].func.main(args[1])
    except AttributeError as e:
        parser.print_help()
        print("Error: No command specified")
        raise


if __name__ == '__main__':
    main()
