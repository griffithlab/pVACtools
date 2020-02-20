import sys
import os
import csv
from lib.prediction_class import *

def split_algorithms(prediction_algorithms):
    class_i_prediction_algorithms = []
    class_ii_prediction_algorithms = []
    for prediction_algorithm in sorted(prediction_algorithms):
        prediction_class = globals()[prediction_algorithm]
        prediction_class_object = prediction_class()
        if isinstance(prediction_class_object, MHCI):
            class_i_prediction_algorithms.append(prediction_algorithm)
        elif isinstance(prediction_class_object, MHCII):
            class_ii_prediction_algorithms.append(prediction_algorithm)
    return (class_i_prediction_algorithms, class_ii_prediction_algorithms)

def split_alleles(alleles):
    class_i_alleles = []
    class_ii_alleles = []
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
    return (class_i_alleles, class_ii_alleles)

def combine_reports(input_files, output_file):
    fieldnames = []
    for input_file in input_files:
        with open(input_file, 'r') as input_file_handle:
            reader = csv.DictReader(input_file_handle, delimiter='\t')
            if len(fieldnames) == 0:
                fieldnames = reader.fieldnames
            else:
                for fieldname in reader.fieldnames:
                    if fieldname not in fieldnames:
                        fieldnames.append(fieldname)

    with open(output_file, 'w') as fout:
        writer = csv.DictWriter(fout, delimiter="\t", restval='NA', fieldnames=fieldnames)
        writer.writeheader()
        for input_file in input_files:
            with open(input_file, 'r') as input_file_handle:
                reader = csv.DictReader(input_file_handle, delimiter='\t')
                for row in reader:
                    writer.writerow(row)
