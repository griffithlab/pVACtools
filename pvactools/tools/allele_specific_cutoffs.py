import sys
import argparse
from pvactools.lib.prediction_class import PredictionClass

def define_parser():
    parser = argparse.ArgumentParser(
        "pvactools allele_specific_cutoffs",
        description="Show the allele specific cutoffs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-a", "--allele",
        help="The allele to use",
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.allele is None:
        PredictionClass.print_all_allele_cutoffs()
    else:
        print(PredictionClass.cutoff_for_allele(args.allele))

if __name__ == "__main__":
    main()
