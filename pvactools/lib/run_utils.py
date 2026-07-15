import sys
import os
import csv
import binascii
import re
from itertools import islice
import argparse
import pandas as pd
from math import ceil

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
            mt_stop = len(mt_epitopes) - i + flanking_length
            wt_stop = len(wt_epitopes) - i + flanking_length
            break
    mutant_subsequence = mt_peptide[start:mt_stop]
    wildtype_subsequence = wt_peptide[start:wt_stop]
    supported_aas = supported_amino_acids()
    if mutant_subsequence[0] not in supported_aas:
        mutant_subsequence = mutant_subsequence[1:]
    if mutant_subsequence[-1] not in supported_aas:
        mutant_subsequence = mutant_subsequence[0:-1]
    if not all([c in supported_aas for c in mutant_subsequence]):
        print("Warning. Mutant sequence contains unsupported amino acid. Skipping entry {}".format(line['index']))
        return
    return mutant_subsequence, wildtype_subsequence

def get_mutated_frameshift_peptide_with_flanking_sequence(wt_peptide, mt_peptide, flanking_length):
    wt_epitopes = determine_neoepitopes(wt_peptide, flanking_length+1)
    mt_epitopes = determine_neoepitopes(mt_peptide, flanking_length+1)
    for start, (wt_epitope, mt_epitope) in enumerate(zip(list(wt_epitopes.values()), list(mt_epitopes.values()))):
        if wt_epitope != mt_epitope:
            break
    mutant_subsequence = mt_peptide[start:]
    wildtype_subsequence = wt_peptide[start:(start + (2 * flanking_length))]
    supported_aas = supported_amino_acids()
    if mutant_subsequence[0] not in supported_aas:
        mutant_subsequence = mutant_subsequence[1:]
    if mutant_subsequence[-1] not in supported_aas:
        mutant_subsequence = mutant_subsequence[0:-1]
    if not all([c in supported_aas for c in mutant_subsequence]):
        print("Warning. Mutant sequence contains unsupported amino acid. Skipping entry {}".format(line['index']))
        return
    return mutant_subsequence, wildtype_subsequence

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

    if tool in ['pvacseq', 'pvacsplice', 'pvacfuse']:
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

def min_match_count(peptide_length):
    return ceil(peptide_length / 2)

def determine_consecutive_matches_from_left(mt_epitope_seq, wt_epitope_seq):
    consecutive_matches = 0
    for a, b in zip(mt_epitope_seq, wt_epitope_seq):
        if a == b:
            consecutive_matches += 1
        else:
            break
    return consecutive_matches

def determine_consecutive_matches_from_right(mt_epitope_seq, wt_epitope_seq):
    consecutive_matches = 0
    for a, b in zip(reversed(mt_epitope_seq), reversed(wt_epitope_seq)):
        if a == b:
            consecutive_matches += 1
        else:
            break
    return consecutive_matches

def determine_total_matches(mt_epitope_seq, wt_epitope_seq):
    matches = 0
    for a, b in zip(mt_epitope_seq, wt_epitope_seq):
        if a == b:
            matches += 1
    return matches

def first_difference(s1, s2):
    # Zip stops at the shorter string length
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        if c1 != c2:
            return i

def valid_tiers(tool):
    if tool == 'pvacseq':
        return ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "PoorTranscript", "LowExpr", "Anchor", "Subclonal", "ProbPos", "Poor", "NoExpr"]
    elif tool == 'pvacfuse':
        return ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "LowReadSupport", "LowExpr", "Anchor", "ProbPos", "Poor"]
    elif tool == 'pvacsplice':
        return ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "PoorTranscript", "LowExpr", "Anchor", "Subclonal", "ProbPos", "Poor", "NoExpr"]
    elif tool == 'pvacbind':
        return ["Pass", "PoorBinder", "PoorImmunogenicity", "PoorPresentation", "RefMatch", "ProbPos", "Poor"]
