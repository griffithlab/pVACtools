import argparse
import os
import sys
from subprocess import run, DEVNULL, STDOUT

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacview run",
        description="Launch pVACview R shiny application",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('pvacseq_dir', help='pVACseq results directory path (e.g. ~/Downloads/pvacseq_run/MHC_Class_I/)')
    parser.add_argument('--r_path', default='R', help='Location of R to be used for launching the app (e.g. /usr/local/bin/R)')
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)
    arguments = ['{}'.format(args.r_path), "-e", "shiny::runApp('{}')".format(args.pvacseq_dir)]
    response = run(arguments, check=True)

if __name__ == '__main__':
    main()
