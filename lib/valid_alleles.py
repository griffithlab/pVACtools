import sys
import argparse
from lib.prediction_class import *

class ValidAlleles:
    def __init__(self, prediction_algorithm):
        self.prediction_algorithm = prediction_algorithm

    def print_valid_alleles(self):
        if self.prediction_algorithm is None:
            print('\n'.join(sorted(PredictionClass.all_valid_allele_names())))
        else:
            prediction_class = globals()[self.prediction_algorithm]
            print("\n".join(sorted(prediction_class().valid_allele_names())))

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser("%s valid_alleles" % tool)
        parser.add_argument(
            "-p", "--prediction-algorithm",
            choices=PredictionClass.prediction_methods(),
            help="The epitope prediction algorithms to use",
        )
        return parser
