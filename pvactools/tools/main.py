import argparse
import sys
import pkg_resources
from pvactools.tools import *

def define_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    #add subcommands
    download_cwls_parser = subparsers.add_parser(
        "download_cwls",
        help="Download pVACtools CWLs for each tool's main pipeline",
        add_help=False
    )
    download_cwls_parser.set_defaults(func=download_cwls)

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
        print(pkg_resources.get_distribution("pvactools").version)
    else:
        try:
            args[0].func.main(args[1])
        except AttributeError as e:
            parser.print_help()
            print("Error: No command specified")
            sys.exit(-1)


if __name__ == '__main__':
    main()
