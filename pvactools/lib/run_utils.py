import sys
import os
import csv
import binascii
import re
from itertools import islice
import argparse

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

def change_permissions_recursive(path, dir_mode, file_mode):
    for root, dirs, files in os.walk(path, topdown=False):
        for dir in [os.path.join(root,d) for d in dirs]:
            os.chmod(dir, dir_mode)
        for file in [os.path.join(root, f) for f in files]:
            os.chmod(file, file_mode)

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

def split_file(reader, lines):
    i = iter(reader)
    piece = list(islice(i, lines))
    while piece:
        yield piece
        piece = list(islice(i, lines))

def construct_index(count, gene, transcript, variant_type, position):
    return '{}.{}.{}.{}.{}'.format(count, gene, transcript, variant_type, position)

def float_range(minimum, maximum):
    """Return function handle of an argument type function for
       ArgumentParser checking a float range: minimum <= arg <= maximum
         minimum - minimum acceptable argument
         maximum - maximum acceptable argument"""

    # Define the function with default arguments
    def float_range_checker(arg):
        """New Type function for argparse - a float within predefined range."""

        try:
            f = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("must be a floating point number")
        if f < minimum or f > maximum:
            raise argparse.ArgumentTypeError("must be in range [" + str(minimum) + " .. " + str(maximum)+"]")
        return f

    # Return function handle to checking function
    return float_range_checker

def pvacsplice_anchors():
    """Return function handle of an argument type function for
       ArgumentParser checking of the pVACsplice anchors
       checking that the specified criteria are in the list of: ['A', 'D', 'NDA', 'DA', 'N']"""

    # Define the function with default arguments
    def pvacsplice_anchors_checker(arg):
        """New Type function for argparse - a comma-separated list with predefined valid values."""

        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in ['A', 'D', 'NDA', 'DA', 'N']:
                raise argparse.ArgumentTypeError("List element must be one of 'A', 'D', 'NDA', 'DA', 'N', not {}".format(argument))
        return arg_list

    # Return function handle to checking function
    return pvacsplice_anchors_checker

def supported_amino_acids():
    return ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

def determine_neoepitopes(sequence, length):
    epitopes = {}
    for i in range(0, len(sequence)-length+1):
        epitopes[i+1] = sequence[i:i+length]
    return epitopes

def get_mutated_peptide_with_flanking_sequence(wt_peptide, mt_peptide, flanking_length):
    wt_epitopes = determine_neoepitopes(wt_peptide, flanking_length+1)
    mt_epitopes = determine_neoepitopes(mt_peptide, flanking_length+1)
    for i in range(1, len(wt_epitopes)):
        wt_epitope = wt_epitopes[i]
        mt_epitope = mt_epitopes[i]
        if wt_epitope != mt_epitope:
            start = i - 1
            break
    for i, (wt_epitope, mt_epitope) in enumerate(zip(reversed(list(wt_epitopes.values())), reversed(list(mt_epitopes.values())))):
        if wt_epitope != mt_epitope:
            stop = len(mt_epitopes) - i + flanking_length
            break
    mutant_subsequence = mt_peptide[start:stop]
    supported_aas = supported_amino_acids()
    if mutant_subsequence[0] not in supported_aas:
        mutant_subsequence = mutant_subsequence[1:]
    if mutant_subsequence[-1] not in supported_aas:
        mutant_subsequence = mutant_subsequence[0:-1]
    if not all([c in supported_aas for c in mutant_subsequence]):
        print("Warning. Mutant sequence contains unsupported amino acid U. Skipping entry {}".format(line['index']))
        return
    return mutant_subsequence

def get_anchor_positions(hla_allele, epitope_length, allele_specific_anchors, anchor_probabilities, anchor_contribution_threshold, mouse_anchor_positions):
        if allele_specific_anchors and epitope_length in anchor_probabilities and hla_allele in anchor_probabilities[epitope_length]:
            probs = anchor_probabilities[epitope_length][hla_allele]
            positions = []
            total_prob = 0
            for (pos, prob) in sorted(probs.items(), key=lambda x: x[1], reverse=True):
                total_prob += float(prob)
                positions.append(int(pos))
                if total_prob > anchor_contribution_threshold:
                    return positions
        elif allele_specific_anchors and epitope_length in mouse_anchor_positions and hla_allele in mouse_anchor_positions[epitope_length]:
            values = mouse_anchor_positions[epitope_length][hla_allele]
            positions = [pos for pos, val in values.items() if val]
            return positions
                
        return [1, 2, epitope_length - 1 , epitope_length]
