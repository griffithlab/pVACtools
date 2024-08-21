import sys
import argparse

from pvactools.lib.prediction_class import *

class ValidAlgorithms:
    def __init__(self, allele, species):
        self.allele = allele
        self.species = species

    def print_valid_algorithms(self):
        if self.allele is None:
            valid_algorithms = []
            if self.species is None:
                valid_algorithms = PredictionClass.prediction_methods()
            else:
                prediction_algorithms = PredictionClass.prediction_methods()
                for algorithm in prediction_algorithms:
                    cls = globals()[algorithm]
                    alleles = cls().valid_allele_names()
                    for allele in alleles:
                        if cls.species_for_allele(allele) == self.species:
                            valid_algorithms.append(algorithm)
                            break
        else:
            PredictionClass.check_alleles_valid([self.allele])
            if (self.species != None and PredictionClass.species_for_allele(self.allele) != self.species):
                raise Exception("Given species does not match given allele.")
                return
            valid_algorithms = []
            prediction_algorithms = PredictionClass.prediction_methods()
            for algorithm in prediction_algorithms:
                cls = globals()[algorithm]
                alleles = cls().valid_allele_names()
                if (self.allele in alleles) \
                      and (PredictionClass.species_for_allele(self.allele) == self.species \
                           or self.species == None):
                    valid_algorithms.append(algorithm)
        print('\n'.join([a for a in valid_algorithms]))

    @classmethod
    def parser(cls, tool="pvacseq"):
        parser = argparse.ArgumentParser(
            "%s valid_algorithms" % tool,
            description="Show a list of algorithms supported given the specified species and/or allele",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            "-a", "--allele",
            help="Show valid algorithms for the selected allele. "
                 + "For a list of available alleles, use: `{} valid_alleles`.".format(tool),
        )
        parser.add_argument(
            "-s", "--species",
            choices=sorted(set(list(PredictionClass.allele_to_species_map().values())), key=str.casefold),
            help="Show valid algorithms for the selected species only",
        )
        return parser
