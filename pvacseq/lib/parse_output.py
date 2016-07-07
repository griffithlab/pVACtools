import argparse
import csv
import re
import operator
import sys
import os
from math import ceil

def protein_identifier_for_label(key_file):
    tsv_reader = csv.reader(key_file, delimiter='\t')
    pattern = re.compile('>')
    key_hash = {}
    for line in tsv_reader:
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

def parse_input_tsv_file(input_tsv_file):
    tsv_reader = csv.DictReader(input_tsv_file, delimiter='\t')
    tsv_entries = {}
    for line in tsv_reader:
        tsv_entries[line['index']] = line
    return tsv_entries

def match_wildtype_and_mutant_entries(netmhc_results, wt_netmhc_results):
    for key, result in netmhc_results.items():
        (wt_netmhc_result_key, mt_position) = key.split('|', 1)
        if result['variant_type'] == 'missense':
            netmhc_results[key]['wt_epitope_seq'] = wt_netmhc_results[wt_netmhc_result_key][mt_position]['wt_epitope_seq']
            netmhc_results[key]['wt_score']       = wt_netmhc_results[wt_netmhc_result_key][mt_position]['wt_score']
        else:
            wt_results        = wt_netmhc_results[wt_netmhc_result_key]
            mt_epitope_seq    = result['mt_epitope_seq']
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

    return netmhc_results

def parse_netmhc_file(input_netmhc_file, tsv_entries, key_file):
    tsv_reader = csv.reader(input_netmhc_file, delimiter='\t')
    protein_identifier_from_label = protein_identifier_for_label(key_file)
    pattern = re.compile('NetMHC|Protein')
    netmhc_results = {}
    wt_netmhc_results = {}
    for line in tsv_reader:
        if len(line) == 0 or pattern.match(line[0]):
            continue
        protein_label = line[0]
        position      = line[1]
        epitope       = line[2]
        score         = line[3]

        if protein_identifier_from_label[protein_label] is not None:
            protein_identifier = protein_identifier_from_label[protein_label]

        (protein_type, tsv_index) = protein_identifier.split('.', 1)
        if protein_type == 'MT':
            tsv_entry = tsv_entries[tsv_index]
            key = "%s|%s" % (tsv_index, position)
            if key not in netmhc_results:
                netmhc_results[key] = {}
            netmhc_results[key]['mt_score']          = int(score)
            netmhc_results[key]['mt_epitope_seq']    = epitope
            netmhc_results[key]['gene_name']         = tsv_entry['gene_name']
            netmhc_results[key]['amino_acid_change'] = tsv_entry['amino_acid_change']
            netmhc_results[key]['variant_type']      = tsv_entry['variant_type']
            netmhc_results[key]['position']          = position
            netmhc_results[key]['tsv_index']         = tsv_index

        if protein_type == 'WT':
            if tsv_index not in wt_netmhc_results:
                wt_netmhc_results[tsv_index] = {}
            wt_netmhc_results[tsv_index][position] = {}
            wt_netmhc_results[tsv_index][position]['wt_score']       = int(score)
            wt_netmhc_results[tsv_index][position]['wt_epitope_seq'] = epitope

    return match_wildtype_and_mutant_entries(netmhc_results, wt_netmhc_results)

def flatten_netmhc_results(netmhc_results):
    #transform the netmhc_results dictionary into a two-dimensional list
    flattened_netmhc_results = list((
        value['gene_name'],
        value['amino_acid_change'],
        value['position'],
        value['mt_score'],
        value['wt_score'],
        value['wt_epitope_seq'],
        value['mt_epitope_seq'],
        value['tsv_index'],
    ) for value in netmhc_results.values())

    return flattened_netmhc_results

def sort_netmhc_results(flattened_netmhc_results):
    sorted_netmhc_results = sorted(
        flattened_netmhc_results,
        key=lambda flattened_netmhc_results: (flattened_netmhc_results[0], flattened_netmhc_results[1], flattened_netmhc_results[3], " ".join(str(item) for item in flattened_netmhc_results))
    )

    return sorted_netmhc_results

def process_input_netmhc_file(input_netmhc_file, tsv_entries, key_file):
    netmhc_results           = parse_netmhc_file(input_netmhc_file, tsv_entries, key_file)
    flattened_netmhc_results = flatten_netmhc_results(netmhc_results)
    sorted_netmhc_results    = sort_netmhc_results(flattened_netmhc_results)

    return sorted_netmhc_results

def output_headers():
    return['Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 'Transcript', 'Ensembl Gene ID', 'Variant Type', 'Mutation', 'Protein Position', 'Gene Name', 'HLA Allele', 'Peptide Length', 'Sub-peptide Position', 'MT score', 'WT score', 'MT epitope seq', 'WT epitope seq', 'Fold Change']

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq parse_output')
    parser.add_argument('input_netmhc_file', type=argparse.FileType('r'), help='Raw output file from Netmhc',)
    parser.add_argument('input_tsv_file', type=argparse.FileType('r'), help='Input list of variants')
    parser.add_argument('key_file', type=argparse.FileType('r'), help='Key file for lookup of FASTA IDs')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='Parsed output file')
    args = parser.parse_args(args_input)

    tsv_writer = csv.DictWriter(args.output_file, delimiter='\t', fieldnames=output_headers())
    tsv_writer.writeheader()

    basename = os.path.basename(args.input_netmhc_file.name)
    (sample, allele, peptide_length, rest) = basename.split(".", 3)

    tsv_entries    = parse_input_tsv_file(args.input_tsv_file)
    netmhc_results = process_input_netmhc_file(args.input_netmhc_file, tsv_entries, args.key_file)
    for gene_name, variant_aa, position, mt_score, wt_score, wt_epitope_seq, mt_epitope_seq, tsv_index in netmhc_results:
        tsv_entry = tsv_entries[tsv_index]
        if mt_epitope_seq != wt_epitope_seq:
            if wt_epitope_seq == 'NA':
                fold_change = 'NA'
            else:
                fold_change = "%.3f" % (wt_score/mt_score)
            tsv_writer.writerow({
                'Chromosome'          : tsv_entry['chromosome_name'],
                'Start'               : tsv_entry['start'],
                'Stop'                : tsv_entry['stop'],
                'Reference'           : tsv_entry['reference'],
                'Variant'             : tsv_entry['variant'],
                'Transcript'          : tsv_entry['transcript_name'],
                'Ensembl Gene ID'     : tsv_entry['ensembl_gene_id'],
                'Variant Type'        : tsv_entry['variant_type'],
                'Mutation'            : variant_aa,
                'Protein Position'    : tsv_entry['protein_position'],
                'Gene Name'           : gene_name,
                'HLA Allele'          : allele,
                'Peptide Length'      : peptide_length,
                'Sub-peptide Position': position,
                'MT score'            : mt_score,
                'WT score'            : wt_score,
                'MT epitope seq'      : mt_epitope_seq,
                'WT epitope seq'      : wt_epitope_seq,
                'Fold Change'         : fold_change,
            })

if __name__ == '__main__':
    main()
