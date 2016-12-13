import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import argparse
import csv
import re
import operator
import os
from math import ceil
from statistics import median
from lib.prediction_class import *
import yaml

csv.field_size_limit(sys.maxsize)

def parse_input_tsv_file(input_tsv_file):
    tsv_reader = csv.DictReader(input_tsv_file, delimiter='\t')
    tsv_entries = {}
    for line in tsv_reader:
        tsv_entries[line['index']] = line
    return tsv_entries

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

def find_mutation_position(wt_epitope_seq, mt_epitope_seq):
    for i,(wt_aa,mt_aa) in enumerate(zip(wt_epitope_seq,mt_epitope_seq)):
        if wt_aa != mt_aa:
            return i+1
    return 0

def match_wildtype_and_mutant_entry_for_missense(result, mt_position, wt_results):
    #The WT epitope at the same position is the match
    match_position = mt_position
    result['wt_epitope_seq']    = wt_results[match_position]['wt_epitope_seq']
    result['wt_scores']         = wt_results[match_position]['wt_scores']
    result['mutation_position'] = find_mutation_position(result['wt_epitope_seq'], result['mt_epitope_seq'])

def match_wildtype_and_mutant_entry_for_frameshift(result, mt_position, wt_results, previous_result):
    #The WT epitope at the same position is the match
    match_position = mt_position
    #Since the MT sequence is longer than the WT sequence, not all MT epitopes have a match
    if match_position not in wt_results:
        result['wt_epitope_seq'] = 'NA'
        result['wt_scores']      = dict.fromkeys(result['mt_scores'].keys(), 'NA')
        if previous_result['mutation_position'] == 'NA':
            result['mutation_position'] = 'NA'
        elif previous_result['mutation_position'] > 0:
            result['mutation_position'] = previous_result['mutation_position'] - 1
        else:
            result['mutation_position'] = 0
        return

    mt_epitope_seq = result['mt_epitope_seq']
    wt_result      = wt_results[match_position]
    wt_epitope_seq = wt_result['wt_epitope_seq']
    if mt_epitope_seq == wt_epitope_seq:
        #The MT epitope does not overlap the frameshift mutation
        result['wt_epitope_seq']    = wt_result['wt_epitope_seq']
        result['wt_scores']         = wt_result['wt_scores']
        result['mutation_position'] = 'NA'
    else:
        #Determine how many consecutive amino acids are the same between the MT epitope and its matching WT epitope
        consecutive_matches = determine_consecutive_matches_from_left(mt_epitope_seq, wt_epitope_seq)
        if consecutive_matches >= min_match_count(int(result['peptide_length'])):
            #The minimum amino acid match count is met
            result['wt_epitope_seq'] = wt_result['wt_epitope_seq']
            result['wt_scores']      = wt_result['wt_scores']
        else:
            #The minimum amino acid match count is not met
            #Even though there is a matching WT epitope there are not enough overlapping amino acids
            #We don't include the matching WT epitope in the output
            result['wt_epitope_seq'] = 'NA'
            result['wt_scores']      = dict.fromkeys(result['mt_scores'].keys(), 'NA')
        mutation_position = find_mutation_position(wt_epitope_seq, mt_epitope_seq)
        if mutation_position == 1 and int(previous_result['mutation_position']) <= 1:
            #The true mutation position is to the left of the current MT eptiope
            mutation_position = 0
        result['mutation_position'] = mutation_position

def match_wildtype_and_mutant_entry_for_inframe_indel(result, mt_position, wt_results, previous_result, iedb_results_for_wt_iedb_result_key):
    #If the previous WT epitope was matched "from the right" we can just use that position to infer the mutation position and match direction
    if previous_result is not None and previous_result['match_direction'] == 'right':
        best_match_position           = previous_result['wt_epitope_position'] + 1
        result['wt_epitope_position'] = best_match_position
        result['match_direction']     = 'right'
        if previous_result['mutation_position'] > 0:
            result['mutation_position'] = previous_result['mutation_position'] - 1
        else:
            result['mutation_position'] = 0

        #We need to ensure that the matched WT eptiope has enough overlapping amino acids with the MT epitope
        best_match_wt_result = wt_results[str(best_match_position)]
        best_match_count     = determine_consecutive_matches_from_right(result['mt_epitope_seq'], best_match_wt_result['wt_epitope_seq'])
        if best_match_count and best_match_count >= min_match_count(int(result['peptide_length'])):
            #The minimum amino acid match count is met
            result['wt_epitope_seq'] = best_match_wt_result['wt_epitope_seq']
            result['wt_scores']      = best_match_wt_result['wt_scores']
        else:
            #The minimum amino acid match count is not met
            #Even though there is a matching WT epitope there are not enough overlapping amino acids
            #We don't include the matching WT epitope in the output
            result['wt_epitope_seq'] = 'NA'
            result['wt_scores']      = dict.fromkeys(result['mt_scores'].keys(), 'NA')

        return

    #In all other cases the WT epitope at the same position is used as the baseline match
    baseline_best_match_position = mt_position

    #For an inframe insertion the MT sequence is longer than the WT sequence
    #In this case not all MT epitopes might have a baseline match
    if baseline_best_match_position not in wt_results:
        result['wt_epitope_seq'] = 'NA'
        result['wt_scores']      = dict.fromkeys(result['mt_scores'].keys(), 'NA')
        #We then infer the mutation position and match direction from the previous MT epitope
        result['match_direction']= previous_result['match_direction']
        if previous_result['mutation_position'] > 0:
            result['mutation_position'] = previous_result['mutation_position'] - 1
        else:
            result['mutation_position'] = 0
        return

    mt_epitope_seq = result['mt_epitope_seq']
    baseline_best_match_wt_result      = wt_results[baseline_best_match_position]
    baseline_best_match_wt_epitope_seq = baseline_best_match_wt_result['wt_epitope_seq']
    #The MT epitope does not overlap the indel mutation
    if baseline_best_match_wt_epitope_seq == mt_epitope_seq:
        result['wt_epitope_seq']      = baseline_best_match_wt_result['wt_epitope_seq']
        result['wt_scores']           = baseline_best_match_wt_result['wt_scores']
        result['wt_epitope_position'] = int(baseline_best_match_position)
        result['mutation_position']   = 'NA'
        result['match_direction']     = 'left'

    #If there is no previous result or the previous WT epitope was matched "from the left" we start by comparing to the baseline match
    if previous_result is None or previous_result['match_direction'] == 'left':
        best_match_count  = determine_consecutive_matches_from_left(mt_epitope_seq, baseline_best_match_wt_epitope_seq)
        #The alternate best match candidate "from the right" is inferred from the baseline best match position and the indel length
        if result['variant_type'] == 'inframe_ins':
            insertion_length              = len(iedb_results_for_wt_iedb_result_key.keys()) - len(wt_results.keys())
            alternate_best_match_position = int(baseline_best_match_position) - insertion_length
        elif result['variant_type'] == 'inframe_del':
            deletion_length                 = len(wt_results.keys()) - len(iedb_results_for_wt_iedb_result_key.keys())
            alternate_best_match_position   = int(baseline_best_match_position) + deletion_length
        if alternate_best_match_position > 0:
            alternate_best_match_wt_result      = wt_results[str(alternate_best_match_position)]
            alternate_best_match_wt_epitope_seq = alternate_best_match_wt_result['wt_epitope_seq']
            consecutive_matches_from_right      = determine_consecutive_matches_from_right(mt_epitope_seq, alternate_best_match_wt_epitope_seq)
            #We then check if the alternate best match epitope has more matching amino acids than the baseline best match epitope
            #If it does, we pick it as the best match
            if consecutive_matches_from_right > best_match_count:
                best_match_count     = consecutive_matches_from_right
                match_direction      = 'right'
                best_match_position  = alternate_best_match_position
                best_match_wt_result = alternate_best_match_wt_result
            else:
                match_direction      = 'left'
                best_match_position  = baseline_best_match_position
                best_match_wt_result = baseline_best_match_wt_result
        else:
            match_direction      = 'left'
            best_match_position  = baseline_best_match_position
            best_match_wt_result = baseline_best_match_wt_result

        #Now that we have found the matching WT epitope we still need to ensure that it has enough overlapping amino acids
        if best_match_count and best_match_count >= min_match_count(int(result['peptide_length'])):
            #The minimum amino acid match count is met
            result['wt_epitope_seq'] = best_match_wt_result['wt_epitope_seq']
            result['wt_scores']      = best_match_wt_result['wt_scores']
        else:
            #The minimum amino acid match count is not met
            #Even though there is a matching WT epitope there are not enough overlapping amino acids
            #We don't include the matching WT epitope in the output
            result['wt_epitope_seq'] = 'NA'
            result['wt_scores']      = dict.fromkeys(result['mt_scores'].keys(), 'NA')

        result['mutation_position']   = find_mutation_position(baseline_best_match_wt_epitope_seq, mt_epitope_seq)
        result['match_direction']     = match_direction
        result['wt_epitope_position'] = best_match_position

def match_wildtype_and_mutant_entries(iedb_results, wt_iedb_results):
    for key in sorted(iedb_results.keys(), key = lambda x: int(x.split('|')[-1])):
        result = iedb_results[key]
        (wt_iedb_result_key, mt_position) = key.split('|', 1)
        previous_mt_position = str(int(mt_position)-1)
        previous_key = '|'.join([wt_iedb_result_key, previous_mt_position])
        if previous_key in iedb_results:
            previous_result = iedb_results[previous_key]
        else:
            previous_result = None
        wt_results = wt_iedb_results[wt_iedb_result_key]
        if result['variant_type'] == 'missense':
            match_wildtype_and_mutant_entry_for_missense(result, mt_position, wt_results)
        elif result['variant_type'] == 'FS':
             match_wildtype_and_mutant_entry_for_frameshift(result, mt_position, wt_results, previous_result)
        elif result['variant_type'] == 'inframe_ins' or result['variant_type'] == 'inframe_del':
            iedb_results_for_wt_iedb_result_key = dict([(key,value) for key, value in iedb_results.items() if key.startswith(wt_iedb_result_key)])
            match_wildtype_and_mutant_entry_for_inframe_indel(result, mt_position, wt_results, previous_result, iedb_results_for_wt_iedb_result_key)

    return iedb_results

def parse_iedb_file(input_iedb_files, tsv_entries, key_file):
    protein_identifiers_from_label = yaml.load(key_file)
    iedb_results = {}
    wt_iedb_results = {}
    for input_iedb_file in input_iedb_files:
        iedb_tsv_reader = csv.DictReader(input_iedb_file, delimiter='\t')
        (sample, method, remainder) = os.path.basename(input_iedb_file.name).split(".", 2)
        for line in iedb_tsv_reader:
            protein_label  = int(line['seq_num'])
            if 'core_peptide' in line:
                position   = str(int(line['start']) - line['peptide'].find(line['core_peptide']))
            else:
                position   = line['start']
            epitope        = line['peptide']
            score          = line['ic50']
            allele         = line['allele']
            peptide_length = len(epitope)

            if protein_identifiers_from_label[protein_label] is not None:
                protein_identifiers = protein_identifiers_from_label[protein_label]

            for protein_identifier in protein_identifiers:
                (protein_type, tsv_index) = protein_identifier.split('.', 1)
                if protein_type == 'MT':
                    tsv_entry = tsv_entries[tsv_index]
                    key = "%s|%s" % (tsv_index, position)
                    if key not in iedb_results:
                        iedb_results[key] = {}
                        iedb_results[key]['mt_scores']         = {}
                        iedb_results[key]['mt_epitope_seq']    = epitope
                        iedb_results[key]['gene_name']         = tsv_entry['gene_name']
                        iedb_results[key]['amino_acid_change'] = tsv_entry['amino_acid_change']
                        iedb_results[key]['variant_type']      = tsv_entry['variant_type']
                        iedb_results[key]['position']          = position
                        iedb_results[key]['tsv_index']         = tsv_index
                        iedb_results[key]['allele']            = allele
                        iedb_results[key]['peptide_length']    = peptide_length
                    iedb_results[key]['mt_scores'][method] = float(score)

                if protein_type == 'WT':
                    if tsv_index not in wt_iedb_results:
                        wt_iedb_results[tsv_index] = {}
                    if position not in wt_iedb_results[tsv_index]:
                        wt_iedb_results[tsv_index][position] = {}
                        wt_iedb_results[tsv_index][position]['wt_scores']     = {}
                    wt_iedb_results[tsv_index][position]['wt_epitope_seq']    = epitope
                    wt_iedb_results[tsv_index][position]['wt_scores'][method] = float(score)

    return match_wildtype_and_mutant_entries(iedb_results, wt_iedb_results)

def add_summary_metrics(iedb_results):
    iedb_results_with_metrics = {}
    for key, value in iedb_results.items():
        mt_scores     = value['mt_scores']
        best_mt_score = sys.maxsize
        for method, score in mt_scores.items():
            if score < best_mt_score:
                best_mt_score        = score
                best_mt_score_method = method
        value['best_mt_score']          = best_mt_score
        value['corresponding_wt_score'] = value['wt_scores'][best_mt_score_method]
        value['best_mt_score_method']   = best_mt_score_method
        value['median_mt_score']        = median(mt_scores.values())
        wt_scores_with_value = [score for score in value['wt_scores'].values() if score != 'NA']
        if not wt_scores_with_value:
            value['median_wt_score']    = 'NA'
        else:
            value['median_wt_score']    = median(wt_scores_with_value)
        iedb_results_with_metrics[key]  = value

    return iedb_results_with_metrics

def pick_top_results(iedb_results, top_score_metric):
    score_at_position = {}
    for key, value in iedb_results.items():
        (tsv_index, position) = key.split('|', 1)
        if tsv_index not in score_at_position.keys():
            score_at_position[tsv_index] = {}
        if top_score_metric == 'median':
            score_at_position[tsv_index][position] = value['median_mt_score']
        elif top_score_metric == 'lowest':
            score_at_position[tsv_index][position] = value['best_mt_score']

    filtered_iedb_results = {}
    for tsv_index, value in score_at_position.items():
        top_score = sys.maxsize
        for position, score in sorted(value.items(), key=lambda x: x[1]):
            top_score_key = "%s|%s" % (tsv_index, position)
            if iedb_results[top_score_key]['wt_epitope_seq'] != iedb_results[top_score_key]['mt_epitope_seq']:
                filtered_iedb_results[top_score_key] = iedb_results[top_score_key]
                break

    return filtered_iedb_results

def flatten_iedb_results(iedb_results):
    #transform the iedb_results dictionary into a two-dimensional list
    flattened_iedb_results = list((
        value['gene_name'],
        value['amino_acid_change'],
        value['position'],
        value['mutation_position'],
        value['mt_scores'],
        value['wt_scores'],
        value['wt_epitope_seq'],
        value['mt_epitope_seq'],
        value['tsv_index'],
        value['allele'],
        value['peptide_length'],
        value['best_mt_score'],
        value['corresponding_wt_score'],
        value['best_mt_score_method'],
        value['median_mt_score'],
        value['median_wt_score'],
    ) for value in iedb_results.values())

    return flattened_iedb_results

def sort_iedb_results(flattened_iedb_results, top_score_metric):
    if top_score_metric == 'median':
        sorted_iedb_results = sorted(
            flattened_iedb_results,
            key=lambda flattened_iedb_results: (
                flattened_iedb_results[0],
                flattened_iedb_results[1],
                flattened_iedb_results[14],
                " ".join(str(item) for item in flattened_iedb_results),
            )
        )
    elif top_score_metric == 'lowest':
        sorted_iedb_results = sorted(
            flattened_iedb_results,
            key=lambda flattened_iedb_results: (
                flattened_iedb_results[0],
                flattened_iedb_results[1],
                flattened_iedb_results[11],
                " ".join(str(item) for item in flattened_iedb_results),
            )
        )

    return sorted_iedb_results

def process_input_iedb_file(input_iedb_files, tsv_entries, key_file, top_result_per_mutation, top_score_metric):
    iedb_results              = parse_iedb_file(input_iedb_files, tsv_entries, key_file)
    iedb_results_with_metrics = add_summary_metrics(iedb_results)
    if top_result_per_mutation == True:
        filtered_iedb_results  = pick_top_results(iedb_results_with_metrics, top_score_metric)
        flattened_iedb_results = flatten_iedb_results(filtered_iedb_results)
    else:
        flattened_iedb_results = flatten_iedb_results(iedb_results_with_metrics)
    sorted_iedb_results       = sort_iedb_results(flattened_iedb_results, top_score_metric)

    return sorted_iedb_results

def base_headers():
    return[
        'Chromosome',
        'Start',
        'Stop',
        'Reference',
        'Variant',
        'Transcript',
        'Ensembl Gene ID',
        'Variant Type',
        'Mutation',
        'Protein Position',
        'Gene Name',
        'HLA Allele',
        'Peptide Length',
        'Sub-peptide Position',
        'Mutation Position',
        'MT Epitope Seq',
        'WT Epitope Seq',
        'Best MT Score Method',
        'Best MT Score',
        'Corresponding WT Score',
        'Corresponding Fold Change',
        'Tumor DNA Depth',
        'Tumor DNA VAF',
        'Tumor RNA Depth',
        'Tumor RNA VAF',
        'Normal Depth',
        'Normal VAF',
        'Gene Expression',
        'Transcript Expression',
        'Median MT Score',
        'Median WT Score',
        'Median Fold Change',
    ]

def output_headers(methods):
    headers = base_headers()
    for method in methods:
        pretty_method = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)
        headers.append("%s WT Score" % pretty_method)
        headers.append("%s MT Score" % pretty_method)

    return headers

def determine_prediction_methods(input_iedb_files):
    methods = set()
    for input_iedb_file in input_iedb_files:
        (sample, method, remainder) = os.path.basename(input_iedb_file.name).split(".", 2)
        methods.add(method)

    return sorted(list(methods))

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq parse_output')
    parser.add_argument('input_iedb_files', type=argparse.FileType('r'), nargs='+', help='Raw output file from IEDB',)
    parser.add_argument('input_tsv_file', type=argparse.FileType('r'), help='Input list of variants')
    parser.add_argument('key_file', type=argparse.FileType('r'), help='Key file for lookup of FASTA IDs')
    parser.add_argument('output_file', help='Parsed output file')
    parser.add_argument('-t', '--top-result-per-mutation', action='store_true', default=False, help='Output top scoring candidate per allele-length per mutation. Default: False')
    parser.add_argument(
        '-m', '--top-score-metric',
        choices=['lowest', 'median'],
        default='median',
        help="The ic50 scoring metric to use when filtering for the top scoring results. " +
        "lowest: Best MT Score - lowest MT ic50 binding score of all chosen prediction methods. " +
        "median: Median MT Score All Methods - median MT ic50 binding score of all chosen prediction methods. " +
        "Default: median"
    )
    args = parser.parse_args(args_input)

    methods = determine_prediction_methods(args.input_iedb_files)
    tmp_output_file = args.output_file + '.tmp'
    tmp_output_filehandle = open(tmp_output_file, 'w')
    tsv_writer = csv.DictWriter(tmp_output_filehandle, delimiter='\t', fieldnames=output_headers(methods))
    tsv_writer.writeheader()

    tsv_entries  = parse_input_tsv_file(args.input_tsv_file)
    iedb_results = process_input_iedb_file(args.input_iedb_files, tsv_entries, args.key_file, args.top_result_per_mutation, args.top_score_metric)
    for (
        gene_name,
        variant_aa,
        position,
        mutation_position,
        mt_scores,
        wt_scores,
        wt_epitope_seq,
        mt_epitope_seq,
        tsv_index, allele,
        peptide_length,
        best_mt_score,
        corresponding_wt_score,
        best_mt_score_method,
        median_mt_score,
        median_wt_score
    ) in iedb_results:
        tsv_entry = tsv_entries[tsv_index]
        if mt_epitope_seq != wt_epitope_seq:
            if wt_epitope_seq == 'NA':
                corresponding_fold_change = 'NA'
            else:
                corresponding_fold_change = "%.3f" % (corresponding_wt_score/best_mt_score)
            if median_wt_score == 'NA':
                median_fold_change = 'NA'
            else:
                median_fold_change = "%.3f" % (median_wt_score/median_mt_score)
            row = {
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
                'Mutation Position'   : mutation_position,
                'MT Epitope Seq'      : mt_epitope_seq,
                'WT Epitope Seq'      : wt_epitope_seq,
                'Best MT Score Method': PredictionClass.prediction_class_name_for_iedb_prediction_method(best_mt_score_method),
                'Best MT Score'       : best_mt_score,
                'Corresponding WT Score'    : corresponding_wt_score,
                'Corresponding Fold Change' : corresponding_fold_change,
                'Median MT Score'     : median_mt_score,
                'Median WT Score'     : median_wt_score,
                'Median Fold Change'  : median_fold_change,
            }
            for method in methods:
                pretty_method = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)
                row["%s WT Score" % pretty_method] = wt_scores[method]
                row["%s MT Score" % pretty_method] = mt_scores[method]
            if 'gene_expression' in tsv_entry:
                row['Gene Expression'] = tsv_entry['gene_expression']
            if 'transcript_expression' in tsv_entry:
                row['Transcript Expression'] = tsv_entry['transcript_expression']
            if 'normal_depth' in tsv_entry:
                row['Normal Depth'] = tsv_entry['normal_depth']
            if 'normal_vaf' in tsv_entry:
                row['Normal VAF'] = tsv_entry['normal_vaf']
            if 'tdna_depth' in tsv_entry:
                row['Tumor DNA Depth'] = tsv_entry['tdna_depth']
            if 'tdna_vaf' in tsv_entry:
                row['Tumor DNA VAF'] = tsv_entry['tdna_vaf']
            if 'trna_depth' in tsv_entry:
                row['Tumor RNA Depth'] = tsv_entry['trna_depth']
            if 'trna_vaf' in tsv_entry:
                row['Tumor RNA VAF'] = tsv_entry['trna_vaf']
            tsv_writer.writerow(row)

    tmp_output_filehandle.close()
    os.replace(tmp_output_file, args.output_file)
    for file_handle in args.input_iedb_files:
        file_handle.close()
    args.input_tsv_file.close()
    args.key_file.close()

if __name__ == '__main__':
    main()
