import sys
import argparse
from lib.prediction_class import *
from lib.valid_alleles import *

def define_parser():
    parser = argparse.ArgumentParser('pvacfuse valid_alleles')
    parser.add_argument(
        "-p", "--prediction-algorithm",
        choices=PredictionClass.prediction_methods(),
        help="The epitope prediction algorithms to use",
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    ValidAlleles(args.prediction_algorithm).print_valid_alleles()

if __name__ == "__main__":
    main()
