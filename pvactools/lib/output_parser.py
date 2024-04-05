from abc import ABCMeta, abstractmethod
import sys
import csv
import re
import operator
import os
from math import ceil, inf
from statistics import median
import yaml

from pvactools.lib.prediction_class import PredictionClass

csv.field_size_limit(sys.maxsize)

class OutputParser(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        self.input_iedb_files        = kwargs['input_iedb_files']
        self.input_tsv_file          = kwargs['input_tsv_file']
        self.key_file                = kwargs['key_file']
        self.output_file             = kwargs['output_file']
        self.sample_name             = kwargs['sample_name']
        self.add_sample_name         = kwargs.get('add_sample_name_column')
        self.flurry_state            = kwargs.get('flurry_state')

    def parse_input_tsv_file(self):
        with open(self.input_tsv_file, 'r') as reader:
            tsv_reader = csv.DictReader(reader, delimiter='\t')
            tsv_entries = {}
            for line in tsv_reader:
                if line['index'] in tsv_entries:
                    sys.exit('Duplicate TSV indexes')
                tsv_entries[line['index']] = line
            return tsv_entries

    def min_match_count(self, peptide_length):
        return ceil(peptide_length / 2)

    def determine_consecutive_matches_from_left(self, mt_epitope_seq, wt_epitope_seq):
        consecutive_matches = 0
        for a, b in zip(mt_epitope_seq, wt_epitope_seq):
            if a == b:
                consecutive_matches += 1
            else:
                break
        return consecutive_matches

    def determine_consecutive_matches_from_right(self, mt_epitope_seq, wt_epitope_seq):
        consecutive_matches = 0
        for a, b in zip(reversed(mt_epitope_seq), reversed(wt_epitope_seq)):
            if a == b:
                consecutive_matches += 1
            else:
                break
        return consecutive_matches

    def determine_total_matches(self, mt_epitope_seq, wt_epitope_seq):
        matches = 0
        for a, b in zip(mt_epitope_seq, wt_epitope_seq):
            if a == b:
                matches += 1
        return matches

    def find_mutation_position(self, wt_epitope_seq, mt_epitope_seq):
        for i,(wt_aa,mt_aa) in enumerate(zip(wt_epitope_seq,mt_epitope_seq)):
            if wt_aa != mt_aa:
                return i+1
        return 0

    def aa_ins_change_len(self, aa_change):
        aac = aa_change.split('/')
        if aac[0][0] == aac[1][0]:
            return len(aac[1]) - len(aac[0])
        else:
            return len(aac[1]) - len(aac[0]) + 1

    def determine_ins_mut_position_from_previous_result(self, previous_result, mt_epitope_seq, result):
        previous_mutation_position = self.position_to_tuple(previous_result['mutation_position'])
        aa_ins_change_len = self.aa_ins_change_len(result['amino_acid_change'])
        if len(previous_mutation_position) == 2:
            if previous_mutation_position[1] == len(mt_epitope_seq)+1 and aa_ins_change_len > previous_mutation_position[1]-previous_mutation_position[0]+1:
                return '{}-{}'.format(previous_mutation_position[0]-1, previous_mutation_position[1])
            elif previous_mutation_position[0] > 1:
                return '{}-{}'.format(previous_mutation_position[0]-1, previous_mutation_position[1]-1)
            else:
                if previous_mutation_position[1] > 2:
                    return '1-{}'.format(previous_mutation_position[1]-1)
                else:
                    return '1' # choose '1' over '1-1' format
        else:
            if previous_mutation_position[0] > 1:
                if aa_ins_change_len > 1:
                    end = (previous_mutation_position[0] - 1) + (aa_ins_change_len - 1)
                    if end > len(mt_epitope_seq):
                        end = len(mt_epitope_seq)
                    return '{}-{}'.format(previous_mutation_position[0]-1, end)
                else:
                    return '{}'.format(previous_mutation_position[0]-1)
            else:
                return '1'

    def find_ins_mut_position(self, wt_epitope_seq, mt_epitope_seq, aa_change, match_direction):
        if mt_epitope_seq == wt_epitope_seq:
            return None
        aal = self.aa_ins_change_len(aa_change)
        mt_start_pos = None
        mt_end_pos = None
        if match_direction == 'left':
            for i,(wt_aa,mt_aa) in enumerate(zip(wt_epitope_seq,mt_epitope_seq)):
                if wt_aa != mt_aa:
                    mt_start_pos = i+1
                    break
            mt_end_pos = mt_start_pos + aal - 1
        elif match_direction == 'right':
            for i,(wt_aa,mt_aa) in enumerate(zip(wt_epitope_seq[::-1],mt_epitope_seq[::-1])):
                if wt_aa != mt_aa:
                    mt_end_pos = i+1
                    break
            mt_start_pos = mt_end_pos - aal + 1

        if mt_end_pos > len(mt_epitope_seq):
            mt_end_pos = len(mt_epitope_seq)
        if mt_start_pos < 1:
            mt_start_pos = 1
        if mt_start_pos == mt_end_pos:
            return (mt_start_pos,)
        else:
            return (mt_start_pos, mt_end_pos)

    def get_percentiles(self, line, method):
        if method.lower() == 'mhcflurry':
            if self.flurry_state == 'both':
                percentiles = {
                    'percentile': line['percentile'],
                    'mhcflurry_presentation_percentile': line['mhcflurry_presentation_percentile'],
                }
            elif self.flurry_state == 'EL_only':
                percentiles = {'mhcflurry_presentation_percentile': line['mhcflurry_presentation_percentile']}
            else:
                percentiles = {'percentile': line['percentile']}
        elif 'percentile' in line:
            percentiles = {'percentile': line['percentile']}
        elif 'percentile_rank' in line:
            percentiles = {'percentile': line['percentile_rank']}
        elif 'rank' in line:
            percentiles = {'percentile': line['rank']}
        else:
            return {'percentile': 'NA'}

        return dict((k, float(v)) if v != 'None' and v is not None and v != "" else (k, 'NA') for k, v in percentiles.items())

    def get_scores(self, line, method):
        if method.lower() == 'mhcflurry':
            if self.flurry_state == 'both':
                return {
                    'ic50': float(line['ic50']),
                    'mhcflurry_processing_score': float(line['mhcflurry_processing_score']),
                    'mhcflurry_presentation_score': float(line['mhcflurry_presentation_score'])
                }
            elif self.flurry_state == 'EL_only':
                return {
                    'mhcflurry_processing_score': float(line['mhcflurry_processing_score']),
                    'mhcflurry_presentation_score': float(line['mhcflurry_presentation_score'])
                }
            else:
                return {'ic50': float(line['ic50'])}
        elif method.lower() == 'deepimmuno':
            return {'score': float(line['immunogenicity'])}
        elif method.lower() == 'bigmhc_el':
            return {'score': float(line['BigMHC_EL'])}
        elif method.lower() == 'bigmhc_im':
            return {'score': float(line['BigMHC_IM'])}
        elif method.lower() == 'netmhcpan_el':
            return {'score': float(line['score'])}
        elif method.lower() == 'netmhciipan_el':
            return {'score': float(line['score'])}
        else:
            return {'ic50': float(line['ic50'])}

    def format_match_na(self, result, metric):
        return {method: {field: 'NA' for field in fields.keys()} for method, fields in result[f'mt_{metric}s'].items()}

    def position_to_tuple(self, mut_pos):
        # '#-#' -> (#, #); '#' -> (#); or 'NA' -> 'NA'
        if mut_pos == 'NA':
            return mut_pos
        elif '-' in mut_pos:
            d_ind = mut_pos.index('-')
            return (int(mut_pos[0:d_ind]), int(mut_pos[d_ind+1:]))
        else:
            return (int(mut_pos),)

    def match_wildtype_and_mutant_entry_for_missense(self, result, mt_position, wt_results, previous_result):
        #The WT epitope at the same position is the match
        match_position = mt_position
        mt_epitope_seq = result['mt_epitope_seq']
        try:
            wt_result      = wt_results[match_position]
        except:
            import pdb
            pdb.set_trace()
        wt_epitope_seq = wt_result['wt_epitope_seq']
        result['wt_epitope_position'] = match_position
        total_matches  = self.determine_total_matches(mt_epitope_seq, wt_epitope_seq)
        if total_matches >= self.min_match_count(int(result['peptide_length'])):
            result['wt_epitope_seq'] = wt_epitope_seq
            result['wt_scores']      = wt_result['wt_scores']
            result['wt_percentiles'] = wt_result['wt_percentiles']
        else:
            result['wt_epitope_seq'] = 'NA'
            result['wt_scores']      = self.format_match_na(result, 'score')
            result['wt_percentiles'] = self.format_match_na(result, 'percentile')

        if mt_epitope_seq == wt_epitope_seq:
            result['mutation_position'] = 'NA'
        else:
            if previous_result:
                previous_mutation_position = self.position_to_tuple(previous_result['mutation_position'])
                if previous_mutation_position == 'NA':
                    result['mutation_position'] = str(self.find_mutation_position(wt_epitope_seq, mt_epitope_seq))
                elif previous_mutation_position[0] > 0:
                    result['mutation_position'] = str(previous_mutation_position[0] - 1)
                else:
                    result['mutation_position'] = '0'
            else:
                result['mutation_position'] = str(self.find_mutation_position(wt_epitope_seq, mt_epitope_seq))

    def match_wildtype_and_mutant_entry_for_frameshift(self, result, mt_position, wt_results, previous_result):
        #vars for later use
        peptide_length = int(result['peptide_length'])
        #The WT epitope at the same position is the match
        match_position = mt_position
        #Since the MT sequence is longer than the WT sequence, not all MT epitopes have a match
        if match_position not in wt_results:
            result['wt_epitope_seq'] = 'NA'
            result['wt_scores']      = self.format_match_na(result, 'score')
            result['wt_percentiles'] = self.format_match_na(result, 'percentile')
            result['wt_epitope_position'] = 'NA'
            previous_mutation_position = self.position_to_tuple(previous_result['mutation_position'])
            if previous_mutation_position == 'NA':
                result['mutation_position'] = 'NA'
            elif previous_mutation_position[0] == 1:
                result['mutation_position'] = '1-{}'.format(peptide_length)
            else:
                result['mutation_position'] = '{}-{}'.format(previous_mutation_position[0]-1, peptide_length)
            return

        mt_epitope_seq = result['mt_epitope_seq']
        wt_result      = wt_results[match_position]
        wt_epitope_seq = wt_result['wt_epitope_seq']
        if mt_epitope_seq == wt_epitope_seq:
            #The MT epitope does not overlap the frameshift mutation
            result['wt_epitope_seq']    = wt_result['wt_epitope_seq']
            result['wt_scores']         = wt_result['wt_scores']
            result['wt_percentiles']    = wt_result['wt_percentiles']
            result['mutation_position'] = 'NA'
            result['wt_epitope_position'] = 'NA'
        else:
            #Determine how many amino acids are the same between the MT epitope and its matching WT epitope
            total_matches = self.determine_total_matches(mt_epitope_seq, wt_epitope_seq)
            if total_matches >= self.min_match_count(peptide_length):
                #The minimum amino acid match count is met
                result['wt_epitope_seq'] = wt_result['wt_epitope_seq']
                result['wt_scores']      = wt_result['wt_scores']
                result['wt_percentiles'] = wt_result['wt_percentiles']
            else:
                #The minimum amino acid match count is not met
                #Even though there is a matching WT epitope there are not enough overlapping amino acids
                #We don't include the matching WT epitope in the output
                result['wt_epitope_seq'] = 'NA'
                result['wt_scores']      = self.format_match_na(result, 'score')
                result['wt_percentiles'] = self.format_match_na(result, 'percentile')
            mutation_position = self.find_mutation_position(wt_epitope_seq, mt_epitope_seq)
            if mutation_position == peptide_length:
                result['mutation_position'] = '{}'.format(mutation_position)
            else:
                result['mutation_position'] = '{}-{}'.format(mutation_position, peptide_length)
            result['wt_epitope_position'] = match_position

    def match_wildtype_and_mutant_entry_for_inframe_indel(self, result, mt_position, wt_results, previous_result, iedb_results_for_wt_iedb_result_key):
        mt_epitope_seq = result['mt_epitope_seq']
        #If the previous WT epitope was matched "from the right" we can just use that position to infer the mutation position and match direction
        if previous_result is not None and previous_result['match_direction'] == 'right':
            best_match_position           = previous_result['wt_epitope_position'] + 1
            result['wt_epitope_position'] = best_match_position
            result['match_direction']     = 'right'
            result['mutation_position']   = self.determine_ins_mut_position_from_previous_result(previous_result, mt_epitope_seq, result)

            #We need to ensure that the matched WT eptiope has enough overlapping amino acids with the MT epitope
            best_match_wt_result = wt_results[str(best_match_position)]
            total_matches        = self.determine_total_matches(result['mt_epitope_seq'], best_match_wt_result['wt_epitope_seq'])
            if total_matches and total_matches >= self.min_match_count(int(result['peptide_length'])):
                #The minimum amino acid match count is met
                result['wt_epitope_seq'] = best_match_wt_result['wt_epitope_seq']
                result['wt_scores']      = best_match_wt_result['wt_scores']
                result['wt_percentiles'] = best_match_wt_result['wt_percentiles']
            else:
                #The minimum amino acid match count is not met
                #Even though there is a matching WT epitope there are not enough overlapping amino acids
                #We don't include the matching WT epitope in the output
                result['wt_epitope_seq'] = 'NA'
                result['wt_scores']      = self.format_match_na(result, 'score')
                result['wt_percentiles'] = self.format_match_na(result, 'percentile')
            return

        #In all other cases the WT epitope at the same position is used as the baseline match
        baseline_best_match_position = mt_position

        #For an inframe insertion the MT sequence is longer than the WT sequence
        #In this case not all MT epitopes might have a baseline match
        if baseline_best_match_position not in wt_results:
            result['wt_epitope_seq'] = 'NA'
            result['wt_scores']      = self.format_match_na(result, 'score')
            result['wt_percentiles'] = self.format_match_na(result, 'percentile')
            #We then infer the mutation position and match direction from the previous MT epitope
            result['match_direction'] = previous_result['match_direction']
            if previous_result['mutation_position'] == 'NA' or previous_result['mutation_position'] == '1':
                result['mutation_position'] = 'NA'
            else:
                result['mutation_position'] = self.determine_ins_mut_position_from_previous_result(previous_result, mt_epitope_seq, result)
            return

        baseline_best_match_wt_result      = wt_results[baseline_best_match_position]
        baseline_best_match_wt_epitope_seq = baseline_best_match_wt_result['wt_epitope_seq']
        #The MT epitope does not overlap the indel mutation
        if baseline_best_match_wt_epitope_seq == mt_epitope_seq:
            result['wt_epitope_seq']      = baseline_best_match_wt_result['wt_epitope_seq']
            result['wt_scores']           = baseline_best_match_wt_result['wt_scores']
            result['wt_percentiles']      = baseline_best_match_wt_result['wt_percentiles']
            result['wt_epitope_position'] = int(baseline_best_match_position)
            result['mutation_position']   = 'NA'
            result['match_direction']     = 'left'
            return

        #If there is no previous result or the previous WT epitope was matched "from the left" we start by comparing to the baseline match
        if previous_result is None or previous_result['match_direction'] == 'left':
            best_match_count  = self.determine_consecutive_matches_from_left(mt_epitope_seq, baseline_best_match_wt_epitope_seq)
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
                consecutive_matches_from_right      = self.determine_consecutive_matches_from_right(mt_epitope_seq, alternate_best_match_wt_epitope_seq)
                #We then check if the alternate best match epitope has more matching amino acids than the baseline best match epitope
                #If it does, we pick it as the best match
                if consecutive_matches_from_right > best_match_count:
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
            total_matches = self.determine_total_matches(mt_epitope_seq, best_match_wt_result['wt_epitope_seq'])
            if total_matches and total_matches >= self.min_match_count(int(result['peptide_length'])):
                #The minimum amino acid match count is met
                result['wt_epitope_seq'] = best_match_wt_result['wt_epitope_seq']
                result['wt_scores']      = best_match_wt_result['wt_scores']
                result['wt_percentiles'] = best_match_wt_result['wt_percentiles']
            else:
                #The minimum amino acid match count is not met
                #Even though there is a matching WT epitope there are not enough overlapping amino acids
                #We don't include the matching WT epitope in the output
                result['wt_epitope_seq'] = 'NA'
                result['wt_scores']      = self.format_match_na(result, 'score')
                result['wt_percentiles'] = self.format_match_na(result, 'percentile')

            if result['variant_type'] == 'inframe_ins':
                mutation_position = self.find_ins_mut_position(baseline_best_match_wt_epitope_seq, mt_epitope_seq, result['amino_acid_change'], match_direction)
                if mutation_position is None:
                    result['mutation_position'] = 'NA'
                else:
                    if previous_result is None:
                        result['mutation_position'] = '{}-{}'.format(mutation_position[0], mutation_position[1]) if len(mutation_position)==2 else '{}'.format(mutation_position[0])
                    else:
                        if previous_result['mutation_position'] == 'NA':
                            result['mutation_position'] = '{}-{}'.format(mutation_position[0], mutation_position[1]) if len(mutation_position)==2 else '{}'.format(mutation_position[0])
                        else:
                            result['mutation_position'] = self.determine_ins_mut_position_from_previous_result(previous_result, mt_epitope_seq, result)
            elif result['variant_type'] == 'inframe_del':
                mutation_position = self.find_mutation_position(baseline_best_match_wt_epitope_seq, mt_epitope_seq)
                result['mutation_position']   = '{}-{}'.format(mutation_position-1, mutation_position)
            result['match_direction']     = match_direction
            result['wt_epitope_position'] = best_match_position

    def match_wildtype_and_mutant_entries(self, iedb_results, wt_iedb_results):
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
                self.match_wildtype_and_mutant_entry_for_missense(result, mt_position, wt_results, previous_result)
            elif result['variant_type'] == 'FS':
                self.match_wildtype_and_mutant_entry_for_frameshift(result, mt_position, wt_results, previous_result)
            elif result['variant_type'] == 'inframe_ins' or result['variant_type'] == 'inframe_del':
                iedb_results_for_wt_iedb_result_key = dict([(key,value) for key, value in iedb_results.items() if key.startswith(wt_iedb_result_key)])
                self.match_wildtype_and_mutant_entry_for_inframe_indel(result, mt_position, wt_results, previous_result, iedb_results_for_wt_iedb_result_key)

        return iedb_results

    @abstractmethod
    def parse_iedb_file(self, tsv_entries):
        pass

    def get_values_for_summary_metrics(self, result, metric, epitope_type):
        mt_values = dict()
        for (method, values) in result['{}_{}s'.format(epitope_type, metric)].items():
            if metric == 'score':
                mt_values[method] = {field: score for field, score in values.items() if field == 'ic50' and score != 'NA'}
            else:
                mt_values[method] = {field: score for field, score in values.items() if score != 'NA'}
            if not mt_values[method]:
                del mt_values[method]
        return mt_values

    def add_summary_metrics(self, iedb_results):
        ## ***** Important *****
        # Will likely be necessary to rename all of the `best_<etc>` format stuff
        # to explicitly note that it is an aggregate for ic50. Might confuse end users
        # if that is not clear
        iedb_results_with_metrics = {}
        for key, result in iedb_results.items():
            for metric in ['score', 'percentile']:
                mt_values = self.get_values_for_summary_metrics(result, metric, 'mt')
                if not mt_values:
                    result['best_mt_{}'.format(metric)]          = 'NA'
                    result['corresponding_wt_{}'.format(metric)] = 'NA'
                    result['best_mt_{}_method'.format(metric)]   = 'NA'
                    result['median_mt_{}'.format(metric)]        = 'NA'
                else:
                    best_mt_value = sys.maxsize
                    for method in sorted(mt_values.keys()):
                        for value in mt_values[method].values():
                            if value < best_mt_value:
                                best_mt_value        = value
                                best_mt_value_method = method
                    result['best_mt_{}'.format(metric)]          = best_mt_value

                    if metric == 'score':
                        corresponding_wt = result['wt_{}s'.format(metric)][best_mt_value_method]['ic50']
                        calculated_median = median([score['ic50'] for score in mt_values.values()])
                    else:
                        corresponding_wts = list(result['wt_{}s'.format(metric)][best_mt_value_method].values())
                        if len(corresponding_wts) > 1:
                            none_na_corresponding_wts = [p for p in corresponding_wts if p != 'NA']
                            if len(none_na_corresponding_wts) == 0:
                                corresponding_wt = 'NA'
                            else:
                                corresponding_wt = min(none_na_corresponding_wts)
                        else:
                            corresponding_wt = corresponding_wts[0]
                        flattend_percentiles = [percentile for method_results in mt_values.values() for percentile in method_results.values()]
                        calculated_median = median(flattend_percentiles)

                    result['corresponding_wt_{}'.format(metric)] = corresponding_wt
                    result['best_mt_{}_method'.format(metric)]   = best_mt_value_method
                    result['median_mt_{}'.format(metric)]        = calculated_median

                # I can probably combine this with code above
                wt_values = self.get_values_for_summary_metrics(result, metric, 'wt')
                if not wt_values:
                    result['median_wt_{}'.format(metric)]    = 'NA'
                else:
                    if metric == 'score':
                        calculated_median = median([score['ic50'] for score in wt_values.values()])
                    else:
                        flattend_percentiles = [percentile for method_results in wt_values.values() for percentile in method_results.values()]
                        calculated_median = median(flattend_percentiles)
                    result['median_wt_{}'.format(metric)]    = calculated_median

                iedb_results_with_metrics[key]  = result

        return iedb_results_with_metrics

    def flatten_iedb_results(self, iedb_results):
        #transform the iedb_results dictionary into a two-dimensional list
        flattened_iedb_results = []
        for value in iedb_results.values():
            row = []
            for key in (
                'gene_name',
                'amino_acid_change',
                'position',
                'mutation_position',
                'mt_scores',
                'wt_scores',
                'mt_percentiles',
                'wt_percentiles',
                'wt_epitope_seq',
                'mt_epitope_seq',
                'tsv_index',
                'allele',
                'peptide_length',
                'best_mt_score',
                'corresponding_wt_score',
                'best_mt_score_method',
                'median_mt_score',
                'median_wt_score',
                'best_mt_percentile',
                'corresponding_wt_percentile',
                'best_mt_percentile_method',
                'median_mt_percentile',
                'median_wt_percentile',
            ):
                if key in value.keys():
                    row.append(value[key])
                else:
                    row.append('NA')
            flattened_iedb_results.append(row)
        return flattened_iedb_results

    def process_input_iedb_file(self, tsv_entries):
        iedb_results = self.parse_iedb_file(tsv_entries)
        iedb_results_with_metrics = self.add_summary_metrics(iedb_results)
        flattened_iedb_results = self.flatten_iedb_results(iedb_results_with_metrics)

        return flattened_iedb_results

    def base_headers(self):
        headers = [
            'Chromosome',
            'Start',
            'Stop',
            'Reference',
            'Variant',
            'Transcript',
            'Transcript Support Level',
            'Transcript Length',
            'Biotype',
            'Ensembl Gene ID',
            'Variant Type',
            'Mutation',
            'Protein Position',
            'Gene Name',
            'HGVSc',
            'HGVSp',
            'HLA Allele',
            'Peptide Length',
            'Sub-peptide Position',
            'Mutation Position',
            'MT Epitope Seq',
            'WT Epitope Seq',
            'Best MT IC50 Score Method',
            'Best MT IC50 Score',
            'Corresponding WT IC50 Score',
            'Corresponding Fold Change',
            'Best MT Percentile Method',
            'Best MT Percentile',
            'Corresponding WT Percentile',
            'Tumor DNA Depth',
            'Tumor DNA VAF',
            'Tumor RNA Depth',
            'Tumor RNA VAF',
            'Normal Depth',
            'Normal VAF',
            'Gene Expression',
            'Transcript Expression',
            'Median MT IC50 Score',
            'Median WT IC50 Score',
            'Median Fold Change',
            'Median MT Percentile',
            'Median WT Percentile',
        ]
        return headers

    def output_headers(self):
        headers = self.base_headers()
        for method in self.prediction_methods():
            if method.lower() == 'mhcflurry':
                if self.flurry_state == 'EL_only':
                    self.flurry_headers(headers)
                    continue
                elif self.flurry_state == 'both':
                    self.flurry_headers(headers)

            pretty_method = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)
            if method in ['BigMHC_EL', 'BigMHC_IM', 'DeepImmuno']:
                headers.append("%s WT Score" % pretty_method)
                headers.append("%s MT Score" % pretty_method)
                continue

            if method in ['netmhcpan_el', 'netmhciipan_el']:
                headers.append("%s WT Score" % pretty_method)
                headers.append("%s MT Score" % pretty_method)
            else:
                headers.append("%s WT IC50 Score" % pretty_method)
                headers.append("%s MT IC50 Score" % pretty_method)
            headers.append("%s WT Percentile" % pretty_method)
            headers.append("%s MT Percentile" % pretty_method)
        if self.add_sample_name:
            headers.append("Sample Name")
        headers.append("Index")

        return headers

    def flurry_headers(self, headers):
        headers.append("MHCflurryEL Processing WT Score")
        headers.append("MHCflurryEL Processing MT Score")
        headers.append("MHCflurryEL Presentation WT Score")
        headers.append("MHCflurryEL Presentation MT Score")
        headers.append("MHCflurryEL Presentation WT Percentile")
        headers.append("MHCflurryEL Presentation MT Percentile")

    def prediction_methods(self):
        methods = set()
        for input_iedb_file in self.input_iedb_files:
            # we remove "sample_name." prefix from filename and then first part before a dot is the method name 
            method = (os.path.basename(input_iedb_file)[len(self.sample_name)+1:]).split('.', 1)[0]
            methods.add(method)

        return sorted(list(methods))

    def add_pretty_row(self, row, entries, prediction_method, pretty_method, suffix):
        if prediction_method in entries:
            # st = score type, s = score
            for st, s in entries[prediction_method].items():
                if st == 'mhcflurry_presentation_score':
                    row['MHCflurryEL Presentation %s' % suffix.replace("IC50 ", "")] = s
                elif st == 'mhcflurry_processing_score':
                    row['MHCflurryEL Processing %s' % suffix.replace("IC50 ", "")] = s
                elif st == 'mhcflurry_presentation_percentile':
                    row['MHCflurryEL Presentation %s' % suffix.replace("IC50 ", "")] = s
                elif pretty_method in ['NetMHCpanEL', 'NetMHCIIpanEL', 'BigMHC_EL', 'BigMHC_IM', 'DeepImmuno']:
                    row['%s %s' % (pretty_method, suffix.replace("IC50 ", ""))] = s
                else:
                    row['%s %s' % (pretty_method, suffix)] = s
        else:
            # This covers cases where the allele/length of the result row is not supported by the prediction algorithm
            row['%s %s' % (pretty_method, suffix)] = 'NA'

    def execute(self):
        tsv_entries = self.parse_input_tsv_file()
        iedb_results = self.process_input_iedb_file(tsv_entries)

        tmp_output_file = self.output_file + '.tmp'
        tmp_output_filehandle = open(tmp_output_file, 'w')
        tsv_writer = csv.DictWriter(tmp_output_filehandle, delimiter='\t', fieldnames=self.output_headers())
        tsv_writer.writeheader()

        for (
            gene_name,
            variant_aa,
            position,
            mutation_position,
            mt_scores,
            wt_scores,
            mt_percentiles,
            wt_percentiles,
            wt_epitope_seq,
            mt_epitope_seq,
            tsv_index,
            allele,
            peptide_length,
            best_mt_score,
            corresponding_wt_score,
            best_mt_score_method,
            median_mt_score,
            median_wt_score,
            best_mt_percentile,
            corresponding_wt_percentile,
            best_mt_percentile_method,
            median_mt_percentile,
            median_wt_percentile,
        ) in iedb_results:
            tsv_entry = tsv_entries[tsv_index]
            if mt_epitope_seq != wt_epitope_seq:
                if corresponding_wt_score == 'NA':
                    corresponding_fold_change = 'NA'
                elif best_mt_score == 0:
                    corresponding_fold_change = inf
                    corresponding_wt_score = round(corresponding_wt_score, 3)
                else:
                    corresponding_fold_change = round((corresponding_wt_score/best_mt_score), 3)
                    corresponding_wt_score = round(corresponding_wt_score, 3)

                if median_wt_score == 'NA':
                    median_fold_change = 'NA'
                elif median_mt_score == 0:
                    median_fold_change = inf
                    median_wt_score = round(median_wt_score, 3)
                else:
                    median_fold_change = round((median_wt_score/median_mt_score), 3)
                    median_wt_score = round(median_wt_score, 3)
                row = {
                    'Chromosome'          : tsv_entry['chromosome_name'],
                    'Start'               : tsv_entry['start'],
                    'Stop'                : tsv_entry['stop'],
                    'Reference'           : tsv_entry['reference'],
                    'Variant'             : tsv_entry['variant'],
                    'Transcript'          : tsv_entry['transcript_name'],
                    'Transcript Support Level': tsv_entry['transcript_support_level'],
                    'Transcript Length'   : tsv_entry['transcript_length'],
                    'Biotype'             : tsv_entry['biotype'],
                    'Ensembl Gene ID'     : tsv_entry['ensembl_gene_id'],
                    'HGVSc'               : tsv_entry['hgvsc'],
                    'HGVSp'               : tsv_entry['hgvsp'],
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
                    'Corresponding WT IC50 Score': corresponding_wt_score,
                    'Corresponding Fold Change' : corresponding_fold_change,
                    'Median WT IC50 Score'     : median_wt_score,
                    'Median Fold Change'  : median_fold_change,
                    'Index'               : tsv_index,
                }
                row['Best MT IC50 Score Method'] = 'NA' if best_mt_score_method == 'NA' else PredictionClass.prediction_class_name_for_iedb_prediction_method(best_mt_score_method)
                row['Best MT Percentile Method'] = 'NA' if best_mt_percentile_method == 'NA' else PredictionClass.prediction_class_name_for_iedb_prediction_method(best_mt_percentile_method)
                row['Best MT IC50 Score'] = 'NA' if best_mt_score == 'NA' else round(best_mt_score, 3)
                row['Best MT Percentile'] = 'NA' if best_mt_percentile == 'NA' else round(best_mt_percentile, 3)
                row['Corresponding WT Percentile'] = 'NA' if corresponding_wt_percentile == 'NA' else round(corresponding_wt_percentile, 3)
                row['Median MT IC50 Score'] = 'NA' if median_mt_score == 'NA' else round(median_mt_score, 3)
                row['Median MT Percentile'] = 'NA' if median_mt_percentile == 'NA' else round(median_mt_percentile, 3)
                row['Median WT Percentile'] = 'NA' if median_wt_percentile == 'NA' else round(median_wt_percentile, 3)
                for method in self.prediction_methods():
                    pretty_method = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)
                    self.add_pretty_row(row, wt_scores, method, pretty_method, 'WT IC50 Score')
                    self.add_pretty_row(row, mt_scores, method, pretty_method, 'MT IC50 Score')
                    if pretty_method not in ['BigMHC_EL', 'BigMHC_IM', 'DeepImmuno']:
                        self.add_pretty_row(row, wt_percentiles, method, pretty_method, 'WT Percentile')
                        self.add_pretty_row(row, mt_percentiles, method, pretty_method, 'MT Percentile')

                for (tsv_key, row_key) in zip(['gene_expression', 'transcript_expression', 'normal_vaf', 'tdna_vaf', 'trna_vaf'], ['Gene Expression', 'Transcript Expression', 'Normal VAF', 'Tumor DNA VAF', 'Tumor RNA VAF']):
                    if tsv_key in tsv_entry:
                        if tsv_entry[tsv_key] == 'NA':
                            row[row_key] = 'NA'
                        else:
                            row[row_key] = round(float(tsv_entry[tsv_key]), 3)
                for (tsv_key, row_key) in zip(['normal_depth', 'tdna_depth', 'trna_depth'], ['Normal Depth', 'Tumor DNA Depth', 'Tumor RNA Depth']):
                    if tsv_key in tsv_entry:
                        row[row_key] = tsv_entry[tsv_key]
                if self.add_sample_name:
                    row['Sample Name'] = self.sample_name
                tsv_writer.writerow(row)

        tmp_output_filehandle.close()
        os.replace(tmp_output_file, self.output_file)

class DefaultOutputParser(OutputParser):
    def parse_iedb_file(self, tsv_entries):
        with open(self.key_file, 'r') as key_file_reader:
            protein_identifiers_from_label = yaml.load(key_file_reader, Loader=yaml.FullLoader)
        iedb_results = {}
        wt_iedb_results = {}
        for input_iedb_file in self.input_iedb_files:
            with open(input_iedb_file, 'r') as reader:
                iedb_tsv_reader = csv.DictReader(reader, delimiter='\t')
                # we remove "sample_name." prefix from filename and then first part before a dot is the method name 
                method = (os.path.basename(input_iedb_file)[len(self.sample_name)+1:]).split('.', 1)[0]
                for line in iedb_tsv_reader:
                    if "Warning: Potential DNA sequence(s)" in line['allele']:
                        continue
                    protein_label  = int(line['seq_num'])
                    if 'core_peptide' in line and int(line['end']) - int(line['start']) == 8:
                        #Start and end refer to the position of the core peptide
                        #Infer the (start) position of the peptide from the positions of the core peptide
                        position   = str(int(line['start']) - line['peptide'].find(line['core_peptide']))
                    else:
                        position   = line['start']
                    percentiles    = self.get_percentiles(line, method)
                    epitope        = line['peptide']
                    scores         = self.get_scores(line, method)
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
                                iedb_results[key]                      = {}
                                iedb_results[key]['mt_scores']         = {}
                                iedb_results[key]['mt_percentiles']    = {}
                                iedb_results[key]['mt_epitope_seq']    = epitope
                                iedb_results[key]['gene_name']         = tsv_entry['gene_name']
                                iedb_results[key]['amino_acid_change'] = tsv_entry['amino_acid_change']
                                iedb_results[key]['variant_type']      = tsv_entry['variant_type']
                                iedb_results[key]['position']          = position
                                iedb_results[key]['tsv_index']         = tsv_index
                                iedb_results[key]['allele']            = allele
                                iedb_results[key]['peptide_length']    = peptide_length
                            iedb_results[key]['mt_scores'][method] = scores
                            iedb_results[key]['mt_percentiles'][method] = percentiles
                        else:
                            if tsv_index not in wt_iedb_results:
                                wt_iedb_results[tsv_index] = {}
                            if position not in wt_iedb_results[tsv_index]:
                                wt_iedb_results[tsv_index][position] = {}
                                wt_iedb_results[tsv_index][position][protein_type.lower() + '_scores'] = {}
                                wt_iedb_results[tsv_index][position][protein_type.lower() + '_percentiles'] = {}
                            wt_iedb_results[tsv_index][position][protein_type.lower() + '_epitope_seq'] = epitope
                            wt_iedb_results[tsv_index][position][protein_type.lower() + '_scores'][method] = scores
                            wt_iedb_results[tsv_index][position][protein_type.lower() + '_percentiles'][method] = percentiles

        return self.match_wildtype_and_mutant_entries(iedb_results, wt_iedb_results)


class UnmatchedSequencesOutputParser(OutputParser):
    def parse_iedb_file(self):
        with open(self.key_file, 'r') as key_file_reader:
            tsv_indices_from_label = yaml.load(key_file_reader, Loader=yaml.FullLoader)
        iedb_results = {}
        for input_iedb_file in self.input_iedb_files:
            with open(input_iedb_file, 'r') as reader:
                iedb_tsv_reader = csv.DictReader(reader, delimiter='\t')
                # we remove "sample_name." prefix from filename and then first part before a dot is the method name 
                method = (os.path.basename(input_iedb_file)[len(self.sample_name)+1:]).split('.', 1)[0]
                for line in iedb_tsv_reader:
                    if "Warning: Potential DNA sequence(s)" in line['allele']:
                        continue
                    protein_label  = int(line['seq_num'])
                    if 'core_peptide' in line and int(line['end']) - int(line['start']) == 8:
                        #Start and end refer to the position of the core peptide
                        #Infer the (start) position of the peptide from the positions of the core peptide
                        position   = str(int(line['start']) - line['peptide'].find(line['core_peptide']))
                    else:
                        position   = line['start']
                    epitope        = line['peptide']
                    scores         = self.get_scores(line, method)
                    allele         = line['allele']
                    peptide_length = len(epitope)
                    percentiles    = self.get_percentiles(line, method)

                    if tsv_indices_from_label[protein_label] is not None:
                        tsv_indices = tsv_indices_from_label[protein_label]

                    for index in tsv_indices:
                        key = '|'.join([index, position])
                        if key not in iedb_results:
                            iedb_results[key]                      = {}
                            iedb_results[key]['mt_scores']         = {}
                            iedb_results[key]['mt_percentiles']    = {}
                            iedb_results[key]['mt_epitope_seq']    = epitope
                            iedb_results[key]['position']          = position
                            iedb_results[key]['tsv_index']         = index
                            iedb_results[key]['allele']            = allele
                        iedb_results[key]['mt_scores'][method] = scores
                        iedb_results[key]['mt_percentiles'][method] = percentiles
        return iedb_results

    def add_summary_metrics(self, iedb_results):
        iedb_results_with_metrics = {}
        for key, result in iedb_results.items():
            for metric in ['score', 'percentile']:
                mt_values = self.get_values_for_summary_metrics(result, metric, 'mt')
                if not mt_values:
                    result['best_mt_{}'.format(metric)]          = 'NA'
                    result['best_mt_{}_method'.format(metric)]   = 'NA'
                    result['median_mt_{}'.format(metric)]        = 'NA'
                else:
                    best_mt_value = sys.maxsize
                    for method in sorted(mt_values.keys()):
                        for value in mt_values[method].values():
                            if value < best_mt_value:
                                best_mt_value        = value
                                best_mt_value_method = method
                    result['best_mt_{}'.format(metric)]          = best_mt_value
                    result['best_mt_{}_method'.format(metric)]   = best_mt_value_method

                    if metric == 'score':
                        calculated_median = median([score['ic50'] for score in mt_values.values()])
                    else:
                        flattend_percentiles = [percentile for method_results in mt_values.values() for percentile in method_results.values()]
                        calculated_median = median(flattend_percentiles)
                    result['median_mt_{}'.format(metric)]        = calculated_median
                iedb_results_with_metrics[key]  = result
        return iedb_results_with_metrics

    def flatten_iedb_results(self, iedb_results):
        #transform the iedb_results dictionary into a two-dimensional list
        flattened_iedb_results = list((
            value['position'],
            value['mt_scores'],
            value['mt_percentiles'],
            value['mt_epitope_seq'],
            value['tsv_index'],
            value['allele'],
            value['best_mt_score'],
            value['best_mt_score_method'],
            value['median_mt_score'],
            value['best_mt_percentile'],
            value['best_mt_percentile_method'],
            value['median_mt_percentile'],
        ) for value in iedb_results.values())
        return flattened_iedb_results

    def process_input_iedb_file(self):
        iedb_results              = self.parse_iedb_file()
        iedb_results_with_metrics = self.add_summary_metrics(iedb_results)
        flattened_iedb_results    = self.flatten_iedb_results(iedb_results_with_metrics)
        return flattened_iedb_results

    def base_headers(self):
        return[
            'Mutation',
            'HLA Allele',
            'Sub-peptide Position',
            'Epitope Seq',
            'Median IC50 Score',
            'Best IC50 Score',
            'Best IC50 Score Method',
            'Median Percentile',
            'Best Percentile',
            'Best Percentile Method',
        ]

    def output_headers(self):
        headers = self.base_headers()
        for method in self.prediction_methods():
            if method.lower() == 'mhcflurry':
                if self.flurry_state == 'EL_only':
                    self.flurry_headers(headers)
                    continue 
                elif self.flurry_state == 'both':
                    self.flurry_headers(headers)
            pretty_method = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)
            if method in ['BigMHC_EL', 'BigMHC_IM', 'DeepImmuno']:
                headers.append("%s Score" % pretty_method)
                continue

            if method in ['netmhcpan_el', 'netmhciipan_el']:
                headers.append("%s Score" % pretty_method)
            else:
                headers.append("%s IC50 Score" % pretty_method)
            headers.append("%s Percentile" % pretty_method)
        if self.add_sample_name:
            headers.append("Sample Name")
        return headers

    def flurry_headers(self, headers):
        headers.append("MHCflurryEL Processing Score")
        headers.append("MHCflurryEL Presentation Score")
        headers.append("MHCflurryEL Presentation Percentile")

    def execute(self):
        tmp_output_file = self.output_file + '.tmp'
        tmp_output_filehandle = open(tmp_output_file, 'w')
        tsv_writer = csv.DictWriter(tmp_output_filehandle, delimiter='\t', fieldnames=self.output_headers())
        tsv_writer.writeheader()

        iedb_results = self.process_input_iedb_file()
        for (
            position,
            mt_scores,
            mt_percentiles,
            mt_epitope_seq,
            tsv_index,
            allele,
            best_mt_score,
            best_mt_score_method,
            median_mt_score,
            best_mt_percentile,
            best_mt_percentile_method,
            median_mt_percentile,
        ) in iedb_results:
            row = {
                'HLA Allele'          : allele,
                'Sub-peptide Position': position,
                'Epitope Seq'         : mt_epitope_seq,
                'Best IC50 Score'     : best_mt_score,
                'Best Percentile'     : best_mt_percentile,
                'Mutation'            : tsv_index,
            }
            row['Best IC50 Score Method'] = 'NA' if best_mt_score_method == 'NA' else PredictionClass.prediction_class_name_for_iedb_prediction_method(best_mt_score_method)
            row['Best Percentile Method'] = 'NA' if best_mt_percentile_method == 'NA' else PredictionClass.prediction_class_name_for_iedb_prediction_method(best_mt_percentile_method)
            row['Median IC50 Score'] = 'NA' if median_mt_score == 'NA' else round(median_mt_score, 3)
            row['Median Percentile'] = 'NA' if median_mt_percentile == 'NA' else round(median_mt_percentile, 3)
            for method in self.prediction_methods():
                pretty_method = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)
                self.add_pretty_row(row, mt_scores, method, pretty_method, 'IC50 Score')
                if pretty_method not in ['BigMHC_EL', 'BigMHC_IM', 'DeepImmuno']:
                    self.add_pretty_row(row, mt_percentiles, method, pretty_method, 'Percentile')
            if self.add_sample_name:
                row['Sample Name'] = self.sample_name
            tsv_writer.writerow(row)

        tmp_output_filehandle.close()
        os.replace(tmp_output_file, self.output_file)

