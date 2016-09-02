import argparse
import sys
from subprocess import call
import os
import pkg_resources
try:
    from . import lib
except SystemError:
    import lib

def main():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()

    #add subcommands
    run_main_program_parser = subparsers.add_parser(
        "run",
        help="Runs the pVAC-Seq pipeline",
        add_help=False
    )
    run_main_program_parser.set_defaults(func=lib.main)

    binding_filter_parser = subparsers.add_parser(
        "binding_filter",
        help="Filters variants processed by IEDB by binding score",
        add_help=False
    )
    binding_filter_parser.set_defaults(func=lib.binding_filter)

    coverage_filter_parser = subparsers.add_parser(
        "coverage_filter",
        help="Filters variants processed by IEDB by coverage, vaf, and gene expression",
        add_help=False
    )
    coverage_filter_parser.set_defaults(func=lib.coverage_filter)

    download_example_data_parser = subparsers.add_parser(
        "download_example_data",
        help="Downloads example input and output files",
        add_help=False
    )
    download_example_data_parser.set_defaults(func=lib.download_example_data)

    install_vep_plugin_parser = subparsers.add_parser(
        "install_vep_plugin",
        help="Installs the Wildtype VEP plugin into your VEP_plugins directory",
        add_help=False
    )
    install_vep_plugin_parser.set_defaults(func=lib.install_vep_plugin)

    valid_alleles_parser = subparsers.add_parser(
        "valid_alleles",
        help="Shows a list of valid allele names",
        add_help=False
    )
    valid_alleles_parser.set_defaults(func=lib.valid_alleles)

    config_files_parser = subparsers.add_parser(
        "config_files",
        help="Documentation for the configuration files",
        add_help=False
    )
    config_files_parser.set_defaults(func=lib.config_files)

    parser.add_argument(
        "-v", "--version",
        action="store_true",
        help="Display the currently installed pvacseq version",
    )

    args = parser.parse_known_args()
    if args[0].version is True:
        print(pkg_resources.get_distribution("pvacseq").version)
    else:
        try:
            args[0].func.main(args[1])
        except AttributeError as e:
            parser.print_help()
            print("Error: No command specified")
            sys.exit(-1)


if __name__ == '__main__':
    main()
