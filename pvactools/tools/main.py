import argparse
import sys
try:
    from importlib.metadata import version
except:
    from importlib_metadata import version
from pvactools.tools import *

def define_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    #add subcommands
    allele_specific_cutoffs_parser = subparsers.add_parser(
        "allele_specific_cutoffs",
        help="Show the allele specific cutoffs.",
        add_help=False,
    )
    allele_specific_cutoffs_parser.set_defaults(func=allele_specific_cutoffs)

    compare_parser = subparsers.add_parser(
        "compare",
        help="Run a comparison between two output results folders",
        add_help=False
        )
    compare_parser.set_defaults(func=compare)

    download_cwls_parser = subparsers.add_parser(
        "download_cwls",
        help="Download pVACtools CWLs for each tool's main pipeline",
        add_help=False
    )
    download_cwls_parser.set_defaults(func=download_cwls)

    download_wdls_parser = subparsers.add_parser(
        "download_wdls",
        help="Download pVACtools WDLs to run the main pVACseq and pVACfuse pipelines",
        add_help=False
        )
    download_wdls_parser.set_defaults(func=download_wdls)

    valid_alleles_parser = subparsers.add_parser(
        "valid_alleles",
        help="Show a list of valid allele names.",
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

    parser.add_argument(
        "-v", "--version",
        action="store_true",
        help="Display the currently installed pvactools version",
    )
    return parser

def main():
    parser = define_parser()
    args = parser.parse_known_args()
    if args[0].version is True:
        print(version('pvactools'))
    else:
        try:
            args[0].func.main(args[1])
        except AttributeError as e:
            parser.print_help()
            print("Error: No command specified")
            sys.exit(-1)


if __name__ == '__main__':
    main()
