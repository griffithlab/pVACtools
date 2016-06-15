import argparse
import csv
import re
import operator
import sys
from math import ceil

import pdb

def protein_identifier_for_label(key_file):
    tsvin = csv.reader(key_file, delimiter='\t')
    pattern = re.compile('>')
    key_hash = {}
    for line in tsvin:
        new_name      = line[0]
        original_name = line[1]
        original_name = pattern.sub('', original_name)
        key_hash[new_name] = original_name

    return key_hash

def min_match_count(peptide_length):
    return ceil(peptide_length / 2)

def determine_consecutive_matches(mt_epitope_seq, wt_epitope_seq):
    consecutive_matches = 0
    left_padding        = 0
    #Count consecutive matches from the beginning of the epitope sequences
    for a, b in zip(mt_epitope_seq, wt_epitope_seq):
        if a == b:
            consecutive_matches += 1
            left_padding        += 1
        else:
            break
    #Count consecutive matches from the end of the epitope sequences
    for a, b in zip(reversed(mt_epitope_seq), reversed(wt_epitope_seq)):
        if a == b:
            consecutive_matches += 1
        else:
            break
    return consecutive_matches, left_padding

def parse_input(input_file, key_file):
    tsvin = csv.reader(input_file, delimiter='\t')
    protein_identifier_from_label = protein_identifier_for_label(key_file)
    pattern = re.compile('NetMHC|Protein')
    netmhc_results = {}
    wt_netmhc_results = {}
    for line in tsvin:
        if len(line) == 0 or pattern.match(line[0]):
            continue
        protein_label = line[0]
        position      = line[1]
        epitope       = line[2]
        score         = line[3]

        if protein_identifier_from_label[protein_label] is not None:
            protein_identifier = protein_identifier_from_label[protein_label]

        (protein_type, protein_name, variant_aa) = protein_identifier.split('.', 2)
        if protein_type == 'MT':
            key = "%s_%s" % (protein_identifier, position)
            (tmp, key) = key.split('.', 1)
            if key not in netmhc_results:
                netmhc_results[key] = {}
            netmhc_results[key]['mt_score']       = int(score)
            netmhc_results[key]['mt_epitope_seq'] = epitope
            netmhc_results[key]['protein_name']   = protein_name
            netmhc_results[key]['variant_aa']     = variant_aa
            netmhc_results[key]['position']       = position

        if protein_type == 'WT':
            key = "%s.%s" % (protein_name, variant_aa)
            if key not in wt_netmhc_results:
                wt_netmhc_results[key] = {}
            wt_netmhc_results[key][position] = {}
            wt_netmhc_results[key][position]['wt_score']       = int(score)
            wt_netmhc_results[key][position]['wt_epitope_seq'] = epitope

    for key, result in netmhc_results.items():
        (wt_netmhc_result_key, mt_position) = key.split('_', 1)
        wt_results = wt_netmhc_results[wt_netmhc_result_key]

        mt_epitope_seq = result['mt_epitope_seq']
        best_match_count  = 0
        best_left_padding = 0
        for wt_position, wt_result in wt_results.items():
            wt_epitope_seq = wt_result['wt_epitope_seq']

            consecutive_matches, left_padding = determine_consecutive_matches(mt_epitope_seq, wt_epitope_seq)
            if consecutive_matches > best_match_count:
                best_match_count    = consecutive_matches
                best_left_padding   = left_padding
                best_match_position = wt_position
            elif consecutive_matches == best_match_count and left_padding > best_left_padding:
                best_left_padding   = left_padding
                best_match_position = wt_position

        netmhc_results[key]['wt_epitope_seq'] = 'NA'
        netmhc_results[key]['wt_score']       = 'NA'
        if best_match_count >= min_match_count(len(wt_result['wt_epitope_seq'])):
            netmhc_results[key]['wt_epitope_seq'] = wt_netmhc_results[wt_netmhc_result_key][best_match_position]['wt_epitope_seq']
            netmhc_results[key]['wt_score']       = wt_netmhc_results[wt_netmhc_result_key][best_match_position]['wt_score']

    #transform the netmhc_results dictionary into a two-dimensional list
    netmhc_result_list = list((value['protein_name'], value['variant_aa'], value['position'], value['mt_score'], value['wt_score'], value['wt_epitope_seq'], value['mt_epitope_seq']) for value in netmhc_results.values())
    #sort the list by protein_name, variant_aa, mt_score, and inverse wt_score
    sorted_netmhc_result_list = sorted(netmhc_result_list, key=lambda netmhc_result_list: (netmhc_result_list[0], netmhc_result_list[1], netmhc_result_list[3]))

    return sorted_netmhc_result_list

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq parse_output')
    parser.add_argument('input_file', type=argparse.FileType('r'), help='Raw output file from Netmhc',)
    parser.add_argument('key_file', type=argparse.FileType('r'), help='Key file for lookup of FASTA IDs')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='Parsed output file')
    args = parser.parse_args(args_input)

    tsvout = csv.writer(args.output_file, delimiter='\t', lineterminator='\n')
    tsvout.writerow(['Gene Name', 'Point Mutation', 'Sub-peptide Position', 'MT score', 'WT score', 'MT epitope seq', 'WT epitope seq', 'Fold change'])

    netmhc_results = parse_input(args.input_file, args.key_file)
    for protein_name, variant_aa, position, mt_score, wt_score, wt_epitope_seq, mt_epitope_seq in netmhc_results:
        if mt_epitope_seq != wt_epitope_seq:
            if wt_epitope_seq == 'NA':
                fold_change = 'NA'
            else:
                fold_change = "%.3f" % (wt_score/mt_score)
            tsvout.writerow([protein_name, variant_aa, position, mt_score, wt_score, mt_epitope_seq, wt_epitope_seq, fold_change])


if __name__ == '__main__':
    main()
