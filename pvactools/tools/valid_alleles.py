import sys
import argparse

from pvactools.lib.prediction_class import *

def define_parser():
    parser = argparse.ArgumentParser(
        "pvactools valid_alleles",
        description="Show a list of valid allele names",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-p", "--prediction-algorithm",
        choices=PredictionClass.prediction_methods(),
        help="Show valid alleles for the selected prediction algorithm only",
    )
    parser.add_argument(
        "-s", "--species",
        choices=sorted(set(list(PredictionClass.allele_to_species_map().values())), key=str.casefold),
        help="Show valid alleles for the selected species only",
        default='human'
    )
    parser.add_argument(
        "-l", "--length",
        type=int,
        choices=[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30],
        help="Show valid alleles for the selected length only",
        default=None
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.prediction_algorithm is None:
        if args.length is None:
            length_filtered_alleles = PredictionClass.all_valid_allele_names()
        else:
            length_filtered_alleles = set()
            for algorithm in PredictionClass.prediction_methods():
                prediction_class = globals()[algorithm]()
                alleles = prediction_class.valid_allele_names()
                length_filtered_alleles.update([a for a in alleles if args.length in prediction_class.valid_lengths_for_allele(a)])
    else:
        prediction_class = globals()[args.prediction_algorithm]()
        alleles = prediction_class.valid_allele_names()
        if args.length is None:
            length_filtered_alleles = alleles
        else:
            length_filtered_alleles = [a for a in alleles if args.length in prediction_class.valid_lengths_for_allele(a)]

    print('\n'.join(sorted([a for a in length_filtered_alleles if PredictionClass.species_for_allele(a) == args.species])))

if __name__ == "__main__":
    main()
