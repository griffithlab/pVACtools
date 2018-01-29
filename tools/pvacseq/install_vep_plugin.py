import argparse
import os
import sys
from shutil import copyfile

def define_parser():
    parser = argparse.ArgumentParser('pvacseq install_vep_plugin')
    parser.add_argument('vep_plugins_path', help='Path to your VEP_plugins directory',)
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    base_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    wildtype_plugin_path = os.path.join(base_dir, 'VEP_plugins', 'Wildtype.pm')
    copyfile(wildtype_plugin_path, os.path.join(args.vep_plugins_path, 'Wildtype.pm'))

if __name__ == '__main__':
    main()
