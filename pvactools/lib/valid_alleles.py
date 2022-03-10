import sys
import argparse

from pvactools.lib.prediction_class import *

class ValidAlleles:
    def __init__(self, prediction_algorithm, species):
        self.prediction_algorithm = prediction_algorithm
        self.species = species

    def print_valid_alleles(self):
        if self.prediction_algorithm is None:
            alleles = PredictionClass.all_valid_allele_names()
        else:
            prediction_class = globals()[self.prediction_algorithm]
            alleles = prediction_class().valid_allele_names()
        print('\n'.join(sorted([a for a in alleles if PredictionClass.species_for_allele(a) == self.species])))

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            "%s valid_alleles" % tool,
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
        return parser
