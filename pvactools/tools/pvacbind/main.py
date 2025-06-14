import argparse
import sys
from subprocess import call
import os
from pvactools.tools.pvacbind import *

def define_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    #add subcommands
    run_main_program_parser = subparsers.add_parser(
        "run",
        help="Run the pVACbind pipeline",
        add_help=False
    )
    run_main_program_parser.set_defaults(func=run)

    binding_filter_parser = subparsers.add_parser(
        "binding_filter",
        help="Filter variants processed by IEDB by binding score",
        add_help=False
    )
    binding_filter_parser.set_defaults(func=binding_filter)


    top_score_filter_parser = subparsers.add_parser(
        "top_score_filter",
        help="Pick the best neoepitope for each variant",
        add_help=False,
    )
    top_score_filter_parser.set_defaults(func=top_score_filter)

    net_chop_parser = subparsers.add_parser(
        "net_chop",
        help="Run NetChop on existing pVACbind output .tsv to predict cleavage sites on the neoepitopes.",
        add_help=False,
    )
    net_chop_parser.set_defaults(func=net_chop)

    netmhc_stab_parser = subparsers.add_parser(
        "netmhc_stab",
        help="Run NetMHCStabPan on existing pVACbind output .tsv to add stability predictions to the neoepitopes.",
        add_help=False,
    )
    netmhc_stab_parser.set_defaults(func=netmhc_stab)

    calculate_reference_proteome_similarity_parser = subparsers.add_parser(
        "calculate_reference_proteome_similarity",
        help="Blast peptides against the reference proteome on existing pVACbind output .tsv.",
        add_help=False
    )
    calculate_reference_proteome_similarity_parser.set_defaults(func=calculate_reference_proteome_similarity)

    generate_aggregated_report_parser = subparsers.add_parser(
        "generate_aggregated_report",
        help="Generate an aggregated report from a pVACbind .all_epitopes.tsv report file.",
        add_help=False
    )
    generate_aggregated_report_parser.set_defaults(func=generate_aggregated_report)

    identify_problematic_amino_acids_parser = subparsers.add_parser(
        "identify_problematic_amino_acids",
        help="Mark problematic amino acid positions in each epitope or filter entries that have problematic amino acids.",
        add_help = False
        )
    identify_problematic_amino_acids_parser.set_defaults(func=identify_problematic_amino_acids)

    download_example_data_parser = subparsers.add_parser(
        "download_example_data",
        help="Download example input and output files",
        add_help=False
    )
    download_example_data_parser.set_defaults(func=download_example_data)

    valid_alleles_parser = subparsers.add_parser(
        "valid_alleles",
        help="Show a list of valid allele names",
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
