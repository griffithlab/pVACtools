import sys
import os
import csv

def prediction_method_to_iedb_lookup_dict():
    return{
        'NetMHCpan' : 'netmhcpan',
        'NetMHC'    : 'ann',
        'SMMPMBEC'  : 'smmpmbec',
        'SMM'       : 'smm',
        'NetMHCcons': 'netmhccons',
        'PickPocket': 'pickpocket',
    }

def prediction_methods():
    prediction_method_lookup_dict = prediction_method_to_iedb_lookup_dict()
    return sorted(prediction_method_lookup_dict.keys())

def iedb_to_prediction_method_lookup_dict():
    prediction_method_lookup_dict = prediction_method_to_iedb_lookup_dict()
    return {v: k for k, v in prediction_method_lookup_dict.items()}

def iedb_prediction_methods():
    prediction_method_lookup_dict = iedb_to_prediction_method_lookup_dict()
    return sorted(prediction_method_lookup_dict.keys())

valid_allele_names_for_method_dict = {}

def parse_iedb_allele_file(method):
    #Ultimately we probably want this method to call out to IEDB but their command is currently broken
    #curl --data "method=ann&species=human" http://tools-api.iedb.org/tools_api/mhci/
    base_dir               = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
    iedb_alleles_dir       = os.path.join(base_dir, 'iedb_alleles')
    iedb_alleles_file_name = os.path.join(iedb_alleles_dir, "%s.tsv" % method)
    alleles = {}
    with open(iedb_alleles_file_name) as iedb_alleles_file:
        tsv_reader = csv.DictReader(iedb_alleles_file, delimiter='\t')
        for row in tsv_reader:
            allele = row['MHC']
            if allele not in alleles.keys():
                alleles[allele] = []
            alleles[allele].append(int(row['PeptideLength']))
    return alleles

def valid_allele_names_for_method(method):
    if method not in valid_allele_names_for_method_dict.keys():
        valid_allele_names_for_method_dict[method] = parse_iedb_allele_file(method)

    return valid_allele_names_for_method_dict[method].keys()

def valid_allele_names():
    valid_allele_names = set()
    for method in iedb_prediction_methods():
        valid_allele_names.update(valid_allele_names_for_method(method))
    return list(valid_allele_names)

def check_alleles_valid(alleles):
    valid_alleles = valid_allele_names()
    for allele in alleles:
        if allele not in valid_alleles:
            sys.exit("Allele %s not valid. Run `pvacseq valid_alleles` for a list of valid allele names." % allele)

def check_allele_valid_for_method(allele, method):
    valid_alleles = valid_allele_names_for_method(method)
    if allele not in valid_alleles:
        sys.exit("Allele %s not valid for method %s. Run `pvacseq valid_alleles %s` for a list of valid allele names." % (allele, method, method))

def valid_lengths_for_allele_and_method(allele, method):
    if method not in valid_allele_names_for_method_dict.keys():
        valid_allele_names_for_method_dict[method] = parse_iedb_allele_file(method)

    return valid_allele_names_for_method_dict[method][allele]

def check_length_valid_for_allele_and_method(length, allele, method):
    valid_lengths = valid_lengths_for_allele_and_method(allele, method)
    if length not in valid_lengths:
        sys.exit("Length %s not valid for allele %s and method %s." % (length, allele, method))
