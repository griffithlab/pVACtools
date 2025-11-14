import sys
import os
import csv
import binascii
import re
from itertools import islice
import argparse
import pandas as pd

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

def aggregate_report_evaluations():
    """Return function handle of an argument type function for
       ArgumentParser checking of the aggregate report evaluation values.
       Valid values are: ['Accept', 'Reject', 'Pending', 'Review']"""

    valid_values = ['Accept', 'Reject', 'Pending', 'Review']

    def aggregate_report_evaluation_checker(arg):
        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in valid_values:
                raise argparse.ArgumentTypeError(
                    "Invalid evaluation '{}'. Valid values are: {}".format(argument, ", ".join(valid_values))
                )
        return arg_list

    return aggregate_report_evaluation_checker

def transcript_prioritization_strategy():
    """Return function handle of an argument type function for
       ArgumentParser checking of the transcript prioritization strategy
       checking that the specified criteria are in the list of: ['canonical', 'mane_select', 'tsl']"""

    # Define the function with default arguments
    def transcript_prioritization_strategy_checker(arg):
        """New Type function for argparse - a comma-separated list with predefined valid values."""

        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in ['canonical', 'mane_select', 'tsl']:
                raise argparse.ArgumentTypeError("List element must be one of 'canonical', 'mane_select', 'tsl', not {}".format(argument))
        return arg_list

    # Return function handle to checking function
    return transcript_prioritization_strategy_checker

def top_score_metric2():
    """Return function handle of an argument type function for
       ArgumentParser checking of the top score metric2
       checking that the specified criteria are in the list of: ['ic50', 'combined_percentile', 'binding_percentile', 'immunogenicity_percentile', 'presentation_percentile']"""

    # Define the function with default arguments
    def top_score_metric2_checker(arg):
        """New Type function for argparse - a comma-separated list with predefined valid values."""

        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in ['ic50', 'combined_percentile', 'binding_percentile', 'immunogenicity_percentile', 'presentation_percentile']:
                raise argparse.ArgumentTypeError("List element must be one of 'ic50', 'combined_percentile', 'binding_percentile', 'immunogenicity_percentile', 'presentation_percentile', not {}".format(argument))
        return arg_list

    # Return function handle to checking function
    return top_score_metric2_checker

def pvacsplice_anchors():
    """Return function handle of an argument type function for
       ArgumentParser checking of the pVACsplice anchors
       checking that the specified criteria are in the list of: ['A', 'D', 'NDA']"""

    # Define the function with default arguments
    def pvacsplice_anchors_checker(arg):
        """New Type function for argparse - a comma-separated list with predefined valid values."""

        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in ['A', 'D', 'NDA']:
                raise argparse.ArgumentTypeError("List element must be one of 'A', 'D', 'NDA', not {}".format(argument))
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
    for start, (wt_epitope, mt_epitope) in enumerate(zip(list(wt_epitopes.values()), list(mt_epitopes.values()))):
        if wt_epitope != mt_epitope:
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
        print("Warning. Mutant sequence contains unsupported amino acid. Skipping entry {}".format(line['index']))
        return
    return mutant_subsequence

def get_mutated_frameshift_peptide_with_flanking_sequence(wt_peptide, mt_peptide, flanking_length):
    wt_epitopes = determine_neoepitopes(wt_peptide, flanking_length+1)
    mt_epitopes = determine_neoepitopes(mt_peptide, flanking_length+1)
    for start, (wt_epitope, mt_epitope) in enumerate(zip(list(wt_epitopes.values()), list(mt_epitopes.values()))):
        if wt_epitope != mt_epitope:
            break
    mutant_subsequence = mt_peptide[start:]
    supported_aas = supported_amino_acids()
    if mutant_subsequence[0] not in supported_aas:
        mutant_subsequence = mutant_subsequence[1:]
    if mutant_subsequence[-1] not in supported_aas:
        mutant_subsequence = mutant_subsequence[0:-1]
    if not all([c in supported_aas for c in mutant_subsequence]):
        print("Warning. Mutant sequence contains unsupported amino acid. Skipping entry {}".format(line['index']))
        return
    return mutant_subsequence

def is_preferred_transcript(mutation, transcript_prioritization_strategy, maximum_transcript_support_level):
    if not isinstance(mutation, pd.Series):
        mutation = pd.Series(mutation)
        if mutation['Canonical'] != 'Not Run':
            mutation['Canonical'] = eval(mutation['Canonical'])
        if mutation['MANE Select'] != 'Not Run':
            mutation['MANE Select'] = eval(mutation['MANE Select'])
    if 'mane_select' in transcript_prioritization_strategy:
        if mutation['MANE Select'] == 'Not Run':
            return True
        elif mutation['MANE Select']:
            return True
    if 'canonical' in transcript_prioritization_strategy:
        if mutation['Canonical'] == 'Not Run':
            return True
        elif mutation['Canonical']:
            return True
    if 'tsl' in transcript_prioritization_strategy:
        col = 'TSL' if 'TSL' in mutation else 'Transcript Support Level'
        if pd.isna(mutation[col]):
            return False
        elif mutation[col] == 'NA':
            return False
        elif mutation[col] == 'Not Supported':
            return True
        elif int(mutation[col]) <= maximum_transcript_support_level:
            return True
    return False

def metrics_to_column(tool, metric1, metric2):
    pretty_metric1 = {
        'median': 'Median',
        'lowest': 'Best'
    }
    pretty_metric2 = {
        'ic50': 'IC50 Score',
        'combined_percentile': 'Percentile',
        'binding_percentile': 'IC50 Percentile',
        'immunogenicity_percentile': 'Immunogenicity Percentile',
        'presentation_percentile': 'Presentation Percentile'
    }

    if tool == 'pvacseq':
        return f"{pretty_metric1[metric1]} MT {pretty_metric2[metric2]}"
    else:
        return f"{pretty_metric1[metric1]} {pretty_metric2[metric2]}"

def metric2_to_aggregate_column(metric2):
    pretty_metric2 = {
        'ic50': 'IC50 MT',
        'combined_percentile': '%ile MT',
        'binding_percentile': 'IC50 %ile MT',
        'immunogenicity_percentile': 'IM %ile MT',
        'presentation_percentile': 'Pres %ile MT'
    }
    return pretty_metric2[metric2]
