import argparse
import sys
import os
import pkg_resources
pvac_dir = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
sys.path.append(pvac_dir)
try:
    from . import lib
except (SystemError, ImportError):
    import lib
from tools.pvacfuse import *

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    #add subcommands
    run_main_program_parser = subparsers.add_parser(
        "run",
        help="Runs the pVACfuse pipeline",
        add_help=False
    )
    run_main_program_parser.set_defaults(func=run)

    binding_filter_parser = subparsers.add_parser(
        "binding_filter",
        help="Filters variants processed by IEDB by binding score",
        add_help=False
    )
    binding_filter_parser.set_defaults(func=binding_filter)

    valid_alleles_parser = subparsers.add_parser(
        "valid_alleles",
        help="Shows a list of valid allele names",
        add_help=False
    )
    valid_alleles_parser.set_defaults(func=valid_alleles)

    args = parser.parse_known_args()
    try:
        args[0].func.main(args[1])
    except AttributeError as e:
        parser.print_help()
        print("Error: No command specified")
        sys.exit(-1)

if __name__ == '__main__':
    main()
