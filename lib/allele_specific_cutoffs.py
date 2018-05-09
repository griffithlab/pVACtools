import sys
import argparse
from lib.prediction_class import *

class AlleleSpecificCutoffs:
    def __init__(self, allele):
        self.allele = allele

    def print_allele_specific_cutoffs(self):
        if self.allele is None:
            PredictionClass.print_all_allele_cutoffs()
        else:
            print(PredictionClass.cutoff_for_allele(self.allele))

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser("%s allele_specific_cutoffs" % tool)
        parser.add_argument(
            "-a", "--allele",
            help="The allele to use",
        )
        return parser
