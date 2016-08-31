import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import argparse
from lib.prediction_class import *

def define_parser():
    parser = argparse.ArgumentParser('pvacseq valid_alleles')
    parser.add_argument(
        "-p", "--prediction-algorithm",
        choices=PredictionClass.prediction_methods(),
        help="The epitope prediction algorithms to use",
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.prediction_algorithm is None:
        print('\n'.join(sorted(PredictionClass.all_valid_allele_names())))
    else:
        prediction_class = globals()[args.prediction_algorithm]
        print("\n".join(sorted(prediction_class().valid_allele_names())))

if __name__ == "__main__":
    main()
