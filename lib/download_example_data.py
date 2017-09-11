import argparse
import os
import sys
from shutil import copytree

def define_parser():
    parser = argparse.ArgumentParser('pvacseq download_example_data')
    parser.add_argument('destination_directory', help='Directory for downloading example data',)
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
    source_directory = os.path.join(base_dir, 'example_data')
    copytree(source_directory, os.path.join(args.destination_directory, 'example_data'))

if __name__ == '__main__':
    main()
