from itertools import product

from pvactools.lib.prediction_class import *

def split_algorithms(prediction_algorithms):
    if 'all' in prediction_algorithms:
        return (sorted(MHCI.prediction_methods()), sorted(MHCII.prediction_methods()))
    class_i_prediction_algorithms = set()
    class_ii_prediction_algorithms = set()
    if 'all_class_i' in prediction_algorithms:
        class_i_prediction_algorithms = set(MHCI.prediction_methods())
        prediction_algorithms.remove('all_class_i')
    if 'all_class_ii' in prediction_algorithms:
        class_ii_prediction_algorithms = set(MHCII.prediction_methods())
        prediction_algorithms.remove('all_class_ii')
    for prediction_algorithm in prediction_algorithms:
        prediction_class = globals()[prediction_algorithm]
        prediction_class_object = prediction_class()
        if isinstance(prediction_class_object, MHCI):
            class_i_prediction_algorithms.add(prediction_algorithm)
        elif isinstance(prediction_class_object, MHCII):
            class_ii_prediction_algorithms.add(prediction_algorithm)
    return (sorted(list(class_i_prediction_algorithms)), sorted(list(class_ii_prediction_algorithms)))

def combine_class_ii_alleles(class_ii_alleles):
    valid_combinations = []

    alpha_types = ["DMA", "DOA", "DPA1", "DQA1", "DQA2", "DRA", "DRA", "DRA", "DRA", "DRA"]
    beta_types = ["DMB", "DOB", "DPB1", "DQB1", "DQB2", "DRB1", "DRB2", "DRB3", "DRB4", "DRB5"]
    for (alpha, beta) in zip(alpha_types, beta_types):
        alpha_alleles = []
        beta_alleles = []
        for allele in class_ii_alleles:
            if allele.startswith(alpha):
                alpha_alleles.append(allele)
            elif allele.startswith(beta):
                beta_alleles.append(allele)
        combinations = ["{}-{}".format(alpha_allele, beta_allele) for (alpha_allele, beta_allele) in list(product(alpha_alleles, beta_alleles))]
        valid_combinations.extend([combination for combination in combinations if combination in MHCII.all_valid_allele_names()])

    return sorted(list(set(class_ii_alleles + valid_combinations)))

def split_alleles(alleles):
    class_i_alleles = []
    class_ii_alleles = []
    species = None
    for allele in sorted(set(alleles)):
        valid = 0
        if allele in MHCI.all_valid_allele_names():
            class_i_alleles.append(allele)
            valid = 1
        if allele in MHCII.all_valid_allele_names():
            class_ii_alleles.append(allele)
            valid = 1
        if not valid:
            print("Allele %s not valid. Skipping." % allele)
        else:
            allele_species = PredictionClass.species_for_allele(allele)
            if species is None:
                species = allele_species
            elif species != allele_species:
                raise Exception("Requested alleles are not from the same species.")
    return (class_i_alleles, class_ii_alleles, species)
