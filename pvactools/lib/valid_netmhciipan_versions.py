import argparse


class ValidNetMHCIIPanVersions:
    valid_versions = ["4.3", "4.2", "4.1 (Default)", "4.0 (Not supported by standalone IEDB)"]

    def __init__(self, list):
        self.list = list

    def print_valid_versions(self):
        if self.list:
            print("Valid NetMHCIIpan and NetMHCIIpanEL Versions")
            print('\n'.join([a for a in self.valid_versions]))

    @classmethod
    def parser(cls, tool="pvacseq"):
        parser = argparse.ArgumentParser(
            "%s valid_netmhcpan_versions" % tool,
            description="Show a list of valid versions of NetMHCIIpan and NetMHCIIpanEL that can be used.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            '-l', '--list',
            help="List the valid NetMHCIIpan and NetMHCIIpanEL versions.",
            default='None',
            action='store_true'
        )
        return parser
