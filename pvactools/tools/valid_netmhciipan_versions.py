import sys
import argparse

VALID_VERSIONS = ["4.3", "4.2", "4.1 (Default)", "4.0 (Not supported by standalone IEDB)"]

def define_parser():
    parser = argparse.ArgumentParser(
        "pvactools valid_netmhcpan_versions",
        description="Show a list of valid versions of NetMHCIIpan and NetMHCIIpanEL that can be used.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    print("Valid NetMHCIIpan and NetMHCIIpanEL Versions")
    print('\n'.join([a for a in VALID_VERSIONS]))

if __name__ == "__main__":
    main()
