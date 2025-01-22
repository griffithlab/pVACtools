import sys

from pvactools.lib.netmhc_pan_version import NetMHCPanVersion

def define_parser():
    return NetMHCPanVersion.parser('pvacseq')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    NetMHCPanVersion(args.list).print_valid_versions()

if __name__ == "__main__":
    main()
