import argparse
import sys
from pvactools.tools.pvacvector import *

def define_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    #add subcommands
    run_main_program_parser = subparsers.add_parser(
        "run",
        help="Runs the pVACvector pipeline",
        add_help=False
    )
    run_main_program_parser.set_defaults(func=run)

    run_main_program_parser = subparsers.add_parser(
        "visualize",
        help="Create a visualization of the vector from a pVACvector result fasta",
        add_help=False
    )
    run_main_program_parser.set_defaults(func=visualize)

    valid_alleles_parser = subparsers.add_parser(
        "valid_alleles",
        help="Shows a list of valid allele names",
        add_help=False
    )
    valid_alleles_parser.set_defaults(func=valid_alleles)

    valid_algorithms_parser = subparsers.add_parser(
        "valid_algorithms",
        help="Show a list of algorithms supported given the specified species and/or allele",
        add_help=False
    )
    valid_algorithms_parser.set_defaults(func=valid_algorithms)

    valid_netmhciipan_versions_parser = subparsers.add_parser(
        "valid_netmhciipan_versions",
        help="Show a list of valid versions of NetMHCIIpan and NetMHCIIpanEL that can be used.",
        add_help=False
    )
    valid_netmhciipan_versions_parser.set_defaults(func=valid_netmhciipan_versions)

    allele_specific_cutoffs_parser = subparsers.add_parser(
        "allele_specific_cutoffs",
        help="Show the allele specific cutoffs",
        add_help=False,
    )
    allele_specific_cutoffs_parser.set_defaults(func=allele_specific_cutoffs)

    download_example_data_parser = subparsers.add_parser(
        "download_example_data",
        help="Downloads example input and output files",
        add_help=False
    )
    download_example_data_parser.set_defaults(func=download_example_data)
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
