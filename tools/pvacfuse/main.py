import argparse
import sys
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

    args = parser.parse_known_args()
    try:
        args[0].func.main(args[1])
    except AttributeError as e:
        parser.print_help()
        print("Error: No command specified")
        sys.exit(-1)

if __name__ == '__main__':
    main()
