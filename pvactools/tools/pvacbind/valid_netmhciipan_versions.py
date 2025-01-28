import sys

from pvactools.lib.valid_netmhciipan_versions import ValidNetMHCIIPanVersions

def define_parser():
    return ValidNetMHCIIPanVersions.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    ValidNetMHCIIPanVersions(args.list).print_valid_versions()

if __name__ == "__main__":
    main()
