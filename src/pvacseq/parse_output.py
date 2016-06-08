import argparse
import csv
import re
import operator
import sys

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

def parse_input(input_file, key_file):
    tsvin = csv.reader(input_file, delimiter='\t')
    protein_identifier_from_label = protein_identifier_for_label(key_file)
    pattern = re.compile('NetMHC|Protein')
    netmhc_results = {}
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
        key = "%s.%s" % (protein_identifier, position)
        (tmp, key) = key.split('.', 1)
        if key not in netmhc_results:
            netmhc_results[key] = {}
        if protein_type == 'WT':
            netmhc_results[key]['wt_score']       = int(score)
            netmhc_results[key]['wt_epitope_seq'] = epitope
            netmhc_results[key]['protein_name']   = protein_name
            netmhc_results[key]['variant_aa']     = variant_aa
            netmhc_results[key]['position']       = position
        elif protein_type == 'MT':
            netmhc_results[key]['mt_score']       = int(score)
            netmhc_results[key]['mt_epitope_seq'] = epitope

    #transform the netmhc_results dictionary into a two-dimensional list
    netmhc_result_list = list((value['protein_name'], value['variant_aa'], value['position'], value['mt_score'], value['wt_score'], value['wt_epitope_seq'], value['mt_epitope_seq']) for value in netmhc_results.values())
    #sort the list by protein_name, variant_aa, mt_score, and inverse wt_score
    sorted_netmhc_result_list = sorted(netmhc_result_list, key=lambda netmhc_result_list: (netmhc_result_list[0], netmhc_result_list[1], netmhc_result_list[3], -netmhc_result_list[4]))

    return sorted_netmhc_result_list

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', type=argparse.FileType('r'), help='Raw output file from Netmhc',)
    parser.add_argument('key_file', type=argparse.FileType('r'), help='Key file for lookup of FASTA IDs')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='Parsed output file')
    args = parser.parse_args(args_input)

    tsvout = csv.writer(args.output_file, delimiter='\t', lineterminator='\n')
    tsvout.writerow(['Gene Name', 'Point Mutation', 'Sub-peptide Position', 'MT score', 'WT score', 'MT epitope seq', 'WT epitope seq', 'Fold change'])

    netmhc_results = parse_input(args.input_file, args.key_file)
    for protein_name, variant_aa, position, mt_score, wt_score, wt_epitope_seq, mt_epitope_seq in netmhc_results:
        if mt_epitope_seq != wt_epitope_seq:
            fold_change = "%.3f" % (wt_score/mt_score)
            tsvout.writerow([protein_name, variant_aa, position, mt_score, wt_score, mt_epitope_seq, wt_epitope_seq, fold_change])


if __name__ == '__main__':
    main()
