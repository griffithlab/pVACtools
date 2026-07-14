import sys
import argparse

from pvactools.lib.prediction_class import *

def define_parser():
    parser = argparse.ArgumentParser(
        "pvactools valid_algorithms",
        description="Show a list of algorithms supported given the specified species and/or allele",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-a", "--allele",
        help="Show valid algorithms for the selected allele. "
             + "For a list of available alleles, use: `pvactools valid_alleles`.",
    )
    parser.add_argument(
        "-s", "--species",
        choices=sorted(set(list(PredictionClass.allele_to_species_map().values())), key=str.casefold),
        help="Show valid algorithms for the selected species only",
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.allele is None:
        valid_algorithms = []
        if args.species is None:
            valid_algorithms = PredictionClass.prediction_methods()
        else:
            prediction_algorithms = PredictionClass.prediction_methods()
            for algorithm in prediction_algorithms:
                cls = globals()[algorithm]
                alleles = cls().valid_allele_names()
                for allele in alleles:
                    if cls.species_for_allele(allele) == args.species:
                        valid_algorithms.append(algorithm)
                        break
    else:
        PredictionClass.check_alleles_valid([args.allele])
        if (args.species != None and PredictionClass.species_for_allele(args.allele) != args.species):
            raise Exception("Given species does not match given allele.")
            return
        valid_algorithms = []
        prediction_algorithms = PredictionClass.prediction_methods()
        for algorithm in prediction_algorithms:
            cls = globals()[algorithm]
            alleles = cls().valid_allele_names()
            if (args.allele in alleles) \
                  and (PredictionClass.species_for_allele(args.allele) == args.species \
                       or args.species == None):
                valid_algorithms.append(algorithm)
    print('\n'.join([a for a in valid_algorithms]))

if __name__ == "__main__":
    main()
