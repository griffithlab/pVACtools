import argparse
import os
import sys
from shutil import copytree

class DownloadExampleData:
    def __init__(self, destination_directory, tool):
        self.destination_directory = destination_directory
        self.tool = tool

    def execute(self):
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        source_directory = os.path.join(base_dir, 'tools', self.tool, 'example_data')
        copytree(source_directory, os.path.join(self.destination_directory, '{}_example_data'.format(self.tool)))

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s download_example_data' % tool,
            description="Download example input and output files",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument('destination_directory', help='Directory for downloading example data',)
        return parser
