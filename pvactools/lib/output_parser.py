from abc import ABCMeta, abstractmethod
import sys
import csv
import re
import operator
import requests
import os
import pandas as pd
import h5py
import numpy as np
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
        self.use_normalized_percentiles = kwargs.get('use_normalized_percentiles')
        self.reference_scores_path   = kwargs.get('reference_scores_path')

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

    def find_first_mutation_position(self, wt_epitope_seq, mt_epitope_seq):
        for i,(wt_aa,mt_aa) in enumerate(zip(wt_epitope_seq,mt_epitope_seq)):
            if wt_aa != mt_aa:
                return i+1

    def find_mutation_positions(self, wt_epitope_seq, mt_epitope_seq):
        mutated_positions = []
        for i,(wt_aa,mt_aa) in enumerate(zip(wt_epitope_seq,mt_epitope_seq)):
            if wt_aa != mt_aa:
                mutated_positions.append(i+1)
        if len(mutated_positions) == 0:
            return "NA"
        else:
            return ", ".join([str(x) for x in mutated_positions])

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

    def transform_empty_percentiles(self,p):
        return float(p) if p != 'None' and p is not None and p != "" else 'NA'

    def calculate_normalized_percentile(self, allele, length, score, method, is_reversed=False):
        if allele is None or length is None or score is None or score == 'NA':
            return 'NA'

        normalized = self._normalize_allele(allele)
        if normalized is None:
            return 'NA'

        allele_file = f"HLA-{normalized}_percentiles.h5"
        file_path = os.path.join(self.reference_scores_path, allele_file)

        if not os.path.exists(file_path):
            url = f"https://raw.githubusercontent.com/griffithlab/pvactools_percentiles_data/main/hdf5/{allele_file}"
            try:
                response = requests.get(url, stream=True)
                response.raise_for_status()
                os.makedirs(self.reference_scores_path, exist_ok=True)
                with open(file_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
            except Exception:
                return 'NA'

        key = f"{method}/{length}mer"
        try:
            with h5py.File(file_path, "r") as f:
                if key not in f:
                    return 'NA'  # algorithm or length not present
                ref_scores = f[key][...]
        except Exception:
            return 'NA'

        if ref_scores.size == 0:
            return 'NA'

        n = len(ref_scores)
        left = np.searchsorted(ref_scores, score, side="left")
        right = np.searchsorted(ref_scores, score, side="right")

        if left == right:
            percentile = left / n * 100
        else:
            percentile = (left + right) / (2 * n) * 100

        if is_reversed:
            return round(100.0 - percentile, 3)
        else:
            return round(percentile, 3)

    def _normalize_allele(self, allele):
        # Convert allele to standard format: e.g., HLA-A*01:01 -> HLA-A_01_01
        if allele is None:
            return None

        raw = allele
        if raw.startswith("HLA-"):
            raw = raw[4:]

        try:
            locus, rest = raw.split("*", 1)
        except ValueError:
            return re.sub(r'[^A-Za-z0-9]', '_', raw)

        digits = re.sub(r'\D', '', rest)
        if len(digits) >= 4:
            g1 = digits[0:2]
            g2 = digits[2:]
            normalized = f"{locus}_{g1}_{g2}"
        else:
            normalized = re.sub(r'[^A-Za-z0-9]', '_', raw)

        return normalized

    def _parse_float_or_na(self, value):
        try:
            return float(value)
        except Exception:
            return 'NA'

    def _extract_percentile(self, line, *keys, fallback='NA', is_reversed=False):
        for key in keys:
            if key in line and line[key] not in (None, '', 'NA'):
                return self.transform_empty_percentiles(line[key])
        return fallback

    def _make_score_entry(
        self,
        line,
        label,
        value_key,
        raw_value,
        method,
        percentile_keys=None,
        percentile_fallback='NA',
        include_percentile=True,
        is_reversed=False
    ):
        """
        Centralized entry builder for all predictors.
        Handles:
        - raw numeric parsing
        - percentile extraction
        - normalized percentile override
        """
        val = self._parse_float_or_na(raw_value)
        entry = {label: {value_key: val}}

        if include_percentile:
            if self.use_normalized_percentiles:
                normalized_input = None if val == 'NA' else val
                percentile = self.calculate_normalized_percentile(
                    line.get('allele'),
                    len(line.get('peptide') or ''),
                    normalized_input,
                    method,
                    is_reversed
                )
            else:
                if percentile_keys:
                    percentile = self._extract_percentile(
                        line, *percentile_keys, fallback=percentile_fallback
                    )
                else:
                    percentile = percentile_fallback

            entry[label]['percentile'] = percentile

        return entry

    def get_scores(self, line, method):
        m = method.lower()

        if m == 'mhcflurry':
            if self.flurry_state == 'both':
                return {
                    **self._make_score_entry(
                        line, 'MHCflurry', 'ic50',
                        line.get('ic50'), method,
                        percentile_keys=['percentile']
                    ),
                    **self._make_score_entry(
                        line, 'MHCflurryEL Processing', 'presentation',
                        line.get('mhcflurry_processing_score'), 'MHCflurry_EL_Processing',
                        percentile_keys=None, percentile_fallback='NA', is_reversed=True
                    ),
                    **self._make_score_entry(
                        line, 'MHCflurryEL Presentation', 'presentation',
                        line.get('mhcflurry_presentation_score'), 'MHCflurry_EL_Presentation',
                        percentile_keys=['mhcflurry_presentation_percentile'], is_reversed=True
                    )
                }

            if self.flurry_state == 'el_only':
                return {
                    **self._make_score_entry(
                        line, 'MHCflurryEL Processing', 'presentation',
                        line.get('mhcflurry_processing_score'), 'MHCflurry_EL_Processing',
                        percentile_keys=None, percentile_fallback='NA', is_reversed=True
                    ),
                    **self._make_score_entry(
                        line, 'MHCflurryEL Presentation', 'presentation',
                        line.get('mhcflurry_presentation_score'), 'MHCflurry_EL_Presentation',
                        percentile_keys=['mhcflurry_presentation_percentile'], is_reversed=True
                    )
                }

            return self._make_score_entry(
                line, 'MHCflurry', 'ic50',
                line.get('ic50'), method,
                percentile_keys=['percentile']
            )

        if m == 'deepimmuno':
            return self._make_score_entry(
                line, 'DeepImmuno', 'immunogenicity',
                line.get('immunogenicity'), method,
                percentile_keys=None, percentile_fallback='NA', is_reversed=True
            )

        if m == 'bigmhc_el':
            return self._make_score_entry(
                line, 'BigMHC_EL', 'presentation',
                line.get('BigMHC_EL'), method,
                percentile_keys=None, percentile_fallback='NA', is_reversed=True
            )

        if m == 'bigmhc_im':
            return self._make_score_entry(
                line, 'BigMHC_IM', 'immunogenicity',
                line.get('BigMHC_IM'), method,
                percentile_keys=None, percentile_fallback='NA', is_reversed=True
            )

        if m == 'netmhcpan_el':
            presentation = line.get('score')

            percentile = self._extract_percentile(
                line, 'percentile_rank', 'rank',
                fallback='NA'
            )

            entry = self._make_score_entry(
                line, 'NetMHCpanEL', 'presentation',
                presentation, 'NetMHCpanEL',
                include_percentile=False,
                is_reversed=True
            )

            if not self.use_normalized_percentiles:
                entry['NetMHCpanEL']['percentile'] = percentile
            else:
                entry['NetMHCpanEL']['percentile'] = self.calculate_normalized_percentile(
                    line.get('allele'),
                    len(line.get('peptide') or ''),
                    self._parse_float_or_na(presentation),
                    'NetMHCpanEL',
                    is_reversed=True
                )

            return entry

        if 'netmhciipan_el' in m:
            presentation = (
                line.get('score') if 'score' in line
                else line.get('ic50') if 'ic50' in line
                else None
            )
            if presentation is None:
                raise Exception("Missing expected columns 'score' or 'ic50' in NetMHCIIpanEL output")

            percentile = self._extract_percentile(
                line, 'percentile_rank', 'rank',
                fallback='NA'
            )

            entry = self._make_score_entry(
                line, 'NetMHCIIpanEL', 'presentation',
                presentation, 'NetMHCIIpanEL',
                include_percentile=False,
                is_reversed=True
            )

            if not self.use_normalized_percentiles:
                entry['NetMHCIIpanEL']['percentile'] = percentile
            else:
                entry['NetMHCIIpanEL']['percentile'] = self.calculate_normalized_percentile(
                    line.get('allele'),
                    len(line.get('peptide') or ''),
                    self._parse_float_or_na(presentation),
                    'NetMHCIIpanEL',
                    is_reversed=True
                )

            return entry

        if m == 'mixmhcpred':
            return self._make_score_entry(
                line, 'MixMHCpred', 'binding_score',
                line.get('score'), method,
                percentile_keys=['percentile'],
                is_reversed=True
            )

        if m == 'prime':
            return self._make_score_entry(
                line, 'PRIME', 'immunogenicity',
                line.get('score'), method,
                percentile_keys=['percentile'],
                is_reversed=True
            )

        pretty_method = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)
        percentile_keys = [
            key for key in ('percentile', 'percentile_rank', 'rank') if key in line
        ] or None

        return self._make_score_entry(
            line, pretty_method, 'ic50',
            line.get('ic50'), pretty_method,
            percentile_keys=percentile_keys
        )

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
            result['mutation_position'] = self.find_mutation_positions(wt_epitope_seq, mt_epitope_seq)
        else:
            result['wt_epitope_seq'] = 'NA'
            result['wt_scores']      = self.format_match_na(result, 'score')
            result['mutation_position'] = 'NA'

    def match_wildtype_and_mutant_entry_for_frameshift(self, result, mt_position, wt_results, previous_result):
        #vars for later use
        peptide_length = int(result['peptide_length'])
        #The WT epitope at the same position is the match
        match_position = mt_position
        #Since the MT sequence is longer than the WT sequence, not all MT epitopes have a match
        if match_position not in wt_results:
            result['wt_epitope_seq'] = 'NA'
            result['wt_scores']      = self.format_match_na(result, 'score')
            result['wt_epitope_position'] = 'NA'
            result['mutation_position'] = 'NA'
            return

        mt_epitope_seq = result['mt_epitope_seq']
        wt_result      = wt_results[match_position]
        wt_epitope_seq = wt_result['wt_epitope_seq']
        if mt_epitope_seq == wt_epitope_seq:
            #The MT epitope does not overlap the frameshift mutation
            result['wt_epitope_seq']    = wt_epitope_seq
            result['wt_scores']         = wt_result['wt_scores']
            result['mutation_position'] = 'NA'
            result['wt_epitope_position'] = 'NA'
        else:
            #Determine how many amino acids are the same between the MT epitope and its matching WT epitope
            total_matches = self.determine_total_matches(mt_epitope_seq, wt_epitope_seq)
            if total_matches >= self.min_match_count(peptide_length):
                #The minimum amino acid match count is met
                result['wt_epitope_seq'] = wt_result['wt_epitope_seq']
                result['wt_scores']      = wt_result['wt_scores']
                result['mutation_position'] = self.find_mutation_positions(wt_epitope_seq, mt_epitope_seq)
            else:
                #The minimum amino acid match count is not met
                #Even though there is a matching WT epitope there are not enough overlapping amino acids
                #We don't include the matching WT epitope in the output
                result['wt_epitope_seq'] = 'NA'
                result['wt_scores']      = self.format_match_na(result, 'score')
                result['wt_epitope_seq'] = 'NA'
            result['wt_epitope_position'] = match_position

    def match_wildtype_and_mutant_entry_for_inframe_indel(self, result, mt_position, wt_results, previous_result, iedb_results_for_wt_iedb_result_key):
        mt_epitope_seq = result['mt_epitope_seq']
        #If the previous WT epitope was matched "from the right" we can just use that position to infer the mutation position and match direction
        if previous_result is not None and previous_result['match_direction'] == 'right':
            best_match_position           = previous_result['wt_epitope_position'] + 1
            result['wt_epitope_position'] = best_match_position
            result['match_direction']     = 'right'

            #We need to ensure that the matched WT eptiope has enough overlapping amino acids with the MT epitope
            best_match_wt_result = wt_results[str(best_match_position)]
            total_matches = self.determine_total_matches(result['mt_epitope_seq'], best_match_wt_result['wt_epitope_seq'])
            if total_matches and total_matches >= self.min_match_count(int(result['peptide_length'])):
                #The minimum amino acid match count is met
                result['wt_epitope_seq'] = best_match_wt_result['wt_epitope_seq']
                result['wt_scores']      = best_match_wt_result['wt_scores']
                result['mutation_position'] = self.find_mutation_positions(mt_epitope_seq, best_match_wt_result['wt_epitope_seq'])
            else:
                #The minimum amino acid match count is not met
                #Even though there is a matching WT epitope there are not enough overlapping amino acids
                #We don't include the matching WT epitope in the output
                result['wt_epitope_seq'] = 'NA'
                result['wt_scores']      = self.format_match_na(result, 'score')
                result['mutation_position'] = 'NA'
            return

        #In all other cases the WT epitope at the same position is used as the baseline match
        baseline_best_match_position = mt_position

        #For an inframe insertion the MT sequence is longer than the WT sequence
        #In this case not all MT epitopes might have a baseline match
        if baseline_best_match_position not in wt_results:
            insertion_length = len(iedb_results_for_wt_iedb_result_key.keys()) - len(wt_results.keys())
            best_match_position = int(baseline_best_match_position) - insertion_length
            best_match_wt_result = wt_results[str(best_match_position)]
            result['match_direction'] = 'right'
            result['wt_epitope_position'] = best_match_position
            total_matches = self.determine_total_matches(mt_epitope_seq, best_match_wt_result['wt_epitope_seq'])
            if total_matches and total_matches >= self.min_match_count(int(result['peptide_length'])):
                #The minimum amino acid match count is met
                result['wt_epitope_seq'] = best_match_wt_result['wt_epitope_seq']
                result['wt_scores']      = best_match_wt_result['wt_scores']
                result['mutation_position'] = self.find_mutation_positions(best_match_wt_result['wt_epitope_seq'], mt_epitope_seq)
            else:
                #The minimum amino acid match count is not met
                #Even though there is a matching WT epitope there are not enough overlapping amino acids
                #We don't include the matching WT epitope in the output
                result['wt_epitope_seq'] = 'NA'
                result['wt_scores']      = self.format_match_na(result, 'score')
                result['mutation_position'] = 'NA'
            return

        baseline_best_match_wt_result      = wt_results[baseline_best_match_position]
        baseline_best_match_wt_epitope_seq = baseline_best_match_wt_result['wt_epitope_seq']
        #The MT epitope does not overlap the indel mutation
        if baseline_best_match_wt_epitope_seq == mt_epitope_seq:
            result['wt_epitope_seq']      = baseline_best_match_wt_result['wt_epitope_seq']
            result['wt_scores']           = baseline_best_match_wt_result['wt_scores']
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
                result['mutation_position'] = self.find_mutation_positions(best_match_wt_result['wt_epitope_seq'], mt_epitope_seq)
            else:
                #The minimum amino acid match count is not met
                #Even though there is a matching WT epitope there are not enough overlapping amino acids
                #We don't include the matching WT epitope in the output
                result['wt_epitope_seq'] = 'NA'
                result['wt_scores']      = self.format_match_na(result, 'score')
                result['mutation_position'] = 'NA'

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
        metric_values = dict()
        if metric in ['ic50', 'percentile']:
            for (method, values) in result['{}_scores'.format(epitope_type)].items():
                metric_values[method] = {field: score for field, score in values.items() if field == metric and score != 'NA'}
                if not metric_values[method]:
                    del metric_values[method]
        else:
            if metric == 'ic50_percentile':
                #also include MixMHCpred percentile in ic50_percentile
                result_subset = {method: scores for method, scores in result['{}_scores'.format(epitope_type)].items() if 'ic50' in scores or 'binding_score' in scores}
            else:
                field = metric.replace("_percentile", "")
                result_subset = {method: scores for method, scores in result['{}_scores'.format(epitope_type)].items() if field in scores}
            for (method, values) in result_subset.items():
                metric_values[method] = {field: score for field, score in values.items() if field == 'percentile' and score != 'NA'}
                if not metric_values[method]:
                    del metric_values[method]
        return metric_values

    def add_summary_metrics(self, iedb_results):
        iedb_results_with_metrics = {}
        for key, result in iedb_results.items():
            for metric in ['ic50', 'ic50_percentile', 'immunogenicity_percentile', 'presentation_percentile', 'percentile']:
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
                                best_mt_value = value
                                best_mt_value_method = method
                    result['best_mt_{}'.format(metric)] = best_mt_value
                    result['best_mt_{}_method'.format(metric)]   = best_mt_value_method

                    if metric == 'ic50':
                        result['corresponding_wt_{}'.format(metric)] = result['wt_scores'][best_mt_value_method]['ic50']
                        result['median_mt_{}'.format(metric)] = median([score['ic50'] for score in mt_values.values()])
                    else:
                        result['corresponding_wt_{}'.format(metric)] = result['wt_scores'][best_mt_value_method]['percentile']
                        result['median_mt_{}'.format(metric)] = median([score['percentile'] for score in mt_values.values()])

                wt_values = self.get_values_for_summary_metrics(result, metric, 'wt')
                if not wt_values:
                    result['median_wt_{}'.format(metric)] = 'NA'
                else:
                    if metric == 'ic50':
                        result['median_wt_{}'.format(metric)] = median([score['ic50'] for score in wt_values.values()])
                    else:
                        result['median_wt_{}'.format(metric)] = median([score['percentile'] for score in wt_values.values()])

                iedb_results_with_metrics[key]  = result

        return iedb_results_with_metrics

    def process_input_iedb_file(self, tsv_entries):
        iedb_results = self.parse_iedb_file(tsv_entries)
        iedb_results_with_metrics = self.add_summary_metrics(iedb_results)
        return iedb_results_with_metrics

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
            'Canonical',
            'MANE Select',
            'Biotype',
            'Transcript CDS Flags',
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
            'Best MT IC50 Percentile Method',
            'Best MT IC50 Percentile',
            'Corresponding WT IC50 Percentile',
            'Best MT Immunogenicity Percentile Method',
            'Best MT Immunogenicity Percentile',
            'Corresponding WT Immunogenicity Percentile',
            'Best MT Presentation Percentile Method',
            'Best MT Presentation Percentile',
            'Corresponding WT Presentation Percentile',
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
            'Median MT IC50 Percentile',
            'Median WT IC50 Percentile',
            'Median MT Immunogenicity Percentile',
            'Median WT Immunogenicity Percentile',
            'Median MT Presentation Percentile',
            'Median WT Presentation Percentile',
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
            if method == 'MixMHCpred':
                headers.append("%s WT Binding Score" % pretty_method)
                headers.append("%s MT Binding Score" % pretty_method)
            elif method in ['BigMHC_EL', 'netmhciipan_el', 'netmhcpan_el']:
                headers.append("%s WT Presentation Score" % pretty_method)
                headers.append("%s MT Presentation Score" % pretty_method)
            elif method in ['BigMHC_IM', 'DeepImmuno', 'PRIME']:
                headers.append("%s WT Immunogenicity Score" % pretty_method)
                headers.append("%s MT Immunogenicity Score" % pretty_method)
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
        headers.append("MHCflurryEL Processing WT Percentile")
        headers.append("MHCflurryEL Processing MT Percentile")
        headers.append("MHCflurryEL Presentation WT Score")
        headers.append("MHCflurryEL Presentation MT Score")
        headers.append("MHCflurryEL Presentation WT Percentile")
        headers.append("MHCflurryEL Presentation MT Percentile")

    def prediction_methods(self):
        methods = set()
        pattern = re.compile(rf"{re.escape(self.sample_name)}\.(\w+(?:-\d+\.\d+)?)")

        for input_iedb_file in self.input_iedb_files:
            filename = os.path.basename(input_iedb_file)
            match = pattern.match(filename)
            method = match.group(1)
            methods.add(method)

        return sorted(list(methods))

    def add_prediction_scores(self, row, mt_scores, wt_scores):
        for method in self.prediction_methods():
            pretty_method = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)
            if pretty_method == 'MHCflurry':
                if self.flurry_state == 'EL_only' or self.flurry_state == 'both':
                    row['MHCflurryEL Processing MT Score'] = self.score_or_na(mt_scores, 'MHCflurryEL Processing', 'presentation')
                    row['MHCflurryEL Processing MT Percentile'] = self.score_or_na(mt_scores, 'MHCflurryEL Processing', 'percentile')
                    row['MHCflurryEL Processing WT Score'] = self.score_or_na(wt_scores, 'MHCflurryEL Processing', 'presentation')
                    row['MHCflurryEL Processing WT Percentile'] = self.score_or_na(wt_scores, 'MHCflurryEL Processing', 'percentile')
                    row['MHCflurryEL Presentation MT Score'] = self.score_or_na(mt_scores, 'MHCflurryEL Presentation', 'presentation')
                    row['MHCflurryEL Presentation MT Percentile'] = self.score_or_na(mt_scores, 'MHCflurryEL Presentation', 'percentile')
                    row['MHCflurryEL Presentation WT Score'] = self.score_or_na(wt_scores, 'MHCflurryEL Presentation', 'presentation')
                    row['MHCflurryEL Presentation WT Percentile'] = self.score_or_na(wt_scores, 'MHCflurryEL Presentation', 'percentile')
                if self.flurry_state in ['both', 'BA_only', None]:
                    row['MHCflurry MT IC50 Score'] = self.score_or_na(mt_scores, 'MHCflurry', 'ic50')
                    row['MHCflurry MT Percentile'] = self.score_or_na(mt_scores, 'MHCflurry', 'percentile')
                    row['MHCflurry WT IC50 Score'] = self.score_or_na(wt_scores, 'MHCflurry', 'ic50')
                    row['MHCflurry WT Percentile'] = self.score_or_na(wt_scores, 'MHCflurry', 'percentile')
            else:
                if pretty_method == 'MixMHCpred':
                    row[f'{pretty_method} MT Binding Score'] = self.score_or_na(mt_scores, pretty_method, 'binding_score')
                    row[f'{pretty_method} WT Binding Score'] = self.score_or_na(wt_scores, pretty_method, 'binding_score')
                elif pretty_method in ['BigMHC_EL', 'NetMHCIIpanEL', 'NetMHCpanEL', 'MixMHCpred']:
                    row[f'{pretty_method} MT Presentation Score'] = self.score_or_na(mt_scores, pretty_method, 'presentation')
                    row[f'{pretty_method} WT Presentation Score'] = self.score_or_na(wt_scores, pretty_method, 'presentation')
                elif pretty_method in ['BigMHC_IM', 'DeepImmuno', 'PRIME']:
                    row[f'{pretty_method} MT Immunogenicity Score'] = self.score_or_na(mt_scores, pretty_method, 'immunogenicity')
                    row[f'{pretty_method} WT Immunogenicity Score'] = self.score_or_na(wt_scores, pretty_method, 'immunogenicity')
                else:
                    row[f'{pretty_method} MT IC50 Score'] = self.score_or_na(mt_scores, pretty_method, 'ic50')
                    row[f'{pretty_method} WT IC50 Score'] = self.score_or_na(wt_scores, pretty_method, 'ic50')
                row[f'{pretty_method} MT Percentile'] = self.score_or_na(mt_scores, pretty_method, 'percentile')
                row[f'{pretty_method} WT Percentile'] = self.score_or_na(wt_scores, pretty_method, 'percentile')
        return row

    def score_or_na(self, all_scores, method, score):
        if method in all_scores:
            return all_scores[method][score]
        else:
            return 'NA'

    def rounded_score_or_na(self, score):
        if score == 'NA':
            return score
        else:
            return round(score, 3)

    def execute(self):
        tsv_entries = self.parse_input_tsv_file()
        iedb_results = self.process_input_iedb_file(tsv_entries)

        tmp_output_file = self.output_file + '.tmp'
        tmp_output_filehandle = open(tmp_output_file, 'w')
        tsv_writer = csv.DictWriter(tmp_output_filehandle, delimiter='\t', fieldnames=self.output_headers())
        tsv_writer.writeheader()

        for result in iedb_results.values():
            tsv_entry = tsv_entries[result['tsv_index']]
            if result['mt_epitope_seq'] != result['wt_epitope_seq']:
                if result['corresponding_wt_ic50'] == 'NA':
                    corresponding_fold_change = 'NA'
                elif result['best_mt_ic50'] == 0:
                    corresponding_fold_change = inf
                else:
                    corresponding_fold_change = round((result['corresponding_wt_ic50']/result['best_mt_ic50']), 3)

                if result['median_wt_ic50'] == 'NA':
                    median_fold_change = 'NA'
                elif result['median_mt_ic50'] == 0:
                    median_fold_change = inf
                else:
                    median_fold_change = round((result['median_wt_ic50']/result['median_mt_ic50']), 3)
                row = {
                    'Chromosome'          : tsv_entry['chromosome_name'],
                    'Start'               : tsv_entry['start'],
                    'Stop'                : tsv_entry['stop'],
                    'Reference'           : tsv_entry['reference'],
                    'Variant'             : tsv_entry['variant'],
                    'Transcript'          : tsv_entry['transcript_name'],
                    'Transcript Support Level': tsv_entry['transcript_support_level'],
                    'Transcript Length'   : tsv_entry['transcript_length'],
                    'Canonical'           : tsv_entry['canonical'],
                    'MANE Select'         : tsv_entry['mane_select'],
                    'Biotype'             : tsv_entry['biotype'],
                    'Transcript CDS Flags': tsv_entry['transcript_cds_flags'],
                    'Ensembl Gene ID'     : tsv_entry['ensembl_gene_id'],
                    'HGVSc'               : tsv_entry['hgvsc'],
                    'HGVSp'               : tsv_entry['hgvsp'],
                    'Variant Type'        : tsv_entry['variant_type'],
                    'Mutation'            : result['amino_acid_change'],
                    'Protein Position'    : tsv_entry['protein_position'],
                    'Gene Name'           : result['gene_name'],
                    'HLA Allele'          : result['allele'],
                    'Peptide Length'      : result['peptide_length'],
                    'Sub-peptide Position': result['position'],
                    'Mutation Position'   : result['mutation_position'] if 'mutation_position' in result else 'NA',
                    'MT Epitope Seq'      : result['mt_epitope_seq'],
                    'WT Epitope Seq'      : result['wt_epitope_seq'],
                    'Index'               : result['tsv_index'],
                    #Median IC50 Score
                    'Median MT IC50 Score': self.rounded_score_or_na(result['median_mt_ic50']),
                    'Median WT IC50 Score': self.rounded_score_or_na(result['median_wt_ic50']),
                    'Median Fold Change': median_fold_change,
                    #Median Percentile
                    'Median MT Percentile': self.rounded_score_or_na(result['median_mt_percentile']),
                    'Median WT Percentile': self.rounded_score_or_na(result['median_wt_percentile']),
                    #Median IC50 Percentile
                    'Median MT IC50 Percentile': self.rounded_score_or_na(result['median_mt_ic50_percentile']),
                    'Median WT IC50 Percentile': self.rounded_score_or_na(result['median_wt_ic50_percentile']),
                    #Median Immunogenicity Percentile
                    'Median MT Immunogenicity Percentile': self.rounded_score_or_na(result['median_mt_immunogenicity_percentile']),
                    'Median WT Immunogenicity Percentile': self.rounded_score_or_na(result['median_wt_immunogenicity_percentile']),
                    #Median Presentation Percentile
                    'Median MT Presentation Percentile': self.rounded_score_or_na(result['median_mt_presentation_percentile']),
                    'Median WT Presentation Percentile': self.rounded_score_or_na(result['median_wt_presentation_percentile']),
                    #Best IC50 Score
                    'Best MT IC50 Score': self.rounded_score_or_na(result['best_mt_ic50']),
                    'Best MT IC50 Score Method': result['best_mt_ic50_method'],
                    'Corresponding WT IC50 Score': self.rounded_score_or_na(result['corresponding_wt_ic50']),
                    'Corresponding Fold Change': corresponding_fold_change,
                    #Best Percentile
                    'Best MT Percentile': self.rounded_score_or_na(result['best_mt_percentile']),
                    'Best MT Percentile Method': result['best_mt_percentile_method'],
                    'Corresponding WT Percentile': self.rounded_score_or_na(result['corresponding_wt_percentile']),
                    #Best IC50 Percentile
                    'Best MT IC50 Percentile': self.rounded_score_or_na(result['best_mt_ic50_percentile']),
                    'Best MT IC50 Percentile Method': result['best_mt_ic50_percentile_method'],
                    'Corresponding WT IC50 Percentile': self.rounded_score_or_na(result['corresponding_wt_ic50_percentile']),
                    #Best Immunogenicity Percentile
                    'Best MT Immunogenicity Percentile': self.rounded_score_or_na(result['best_mt_immunogenicity_percentile']),
                    'Best MT Immunogenicity Percentile Method': result['best_mt_immunogenicity_percentile_method'],
                    'Corresponding WT Immunogenicity Percentile': self.rounded_score_or_na(result['corresponding_wt_immunogenicity_percentile']),
                    #Best Presentation Percentile
                    'Best MT Presentation Percentile': self.rounded_score_or_na(result['best_mt_presentation_percentile']),
                    'Best MT Presentation Percentile Method': result['best_mt_presentation_percentile_method'],
                    'Corresponding WT Presentation Percentile': self.rounded_score_or_na(result['corresponding_wt_presentation_percentile']),
                }
                row = self.add_prediction_scores(row, result['mt_scores'], result['wt_scores'])

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
                filename = os.path.basename(input_iedb_file)

                pattern = re.compile(rf"{re.escape(self.sample_name)}\.(\w+(?:-\d+\.\d+)?)")
                match = pattern.match(filename)
                method = match.group(1)

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
                                iedb_results[key]['mt_epitope_seq']    = epitope
                                iedb_results[key]['gene_name']         = tsv_entry['gene_name']
                                iedb_results[key]['amino_acid_change'] = tsv_entry['amino_acid_change']
                                iedb_results[key]['variant_type']      = tsv_entry['variant_type']
                                iedb_results[key]['position']          = position
                                iedb_results[key]['tsv_index']         = tsv_index
                                iedb_results[key]['allele']            = allele
                                iedb_results[key]['peptide_length']    = peptide_length
                            iedb_results[key]['mt_scores'].update(scores)
                        else:
                            if tsv_index not in wt_iedb_results:
                                wt_iedb_results[tsv_index] = {}
                            if position not in wt_iedb_results[tsv_index]:
                                wt_iedb_results[tsv_index][position] = {}
                                wt_iedb_results[tsv_index][position]['wt_scores'] = {}
                            wt_iedb_results[tsv_index][position]['wt_epitope_seq'] = epitope
                            wt_iedb_results[tsv_index][position]['wt_scores'].update(scores)

        return self.match_wildtype_and_mutant_entries(iedb_results, wt_iedb_results)


class UnmatchedSequencesOutputParser(OutputParser):
    def parse_iedb_file(self):
        with open(self.key_file, 'r') as key_file_reader:
            tsv_indices_from_label = yaml.load(key_file_reader, Loader=yaml.FullLoader)
        iedb_results = {}
        for input_iedb_file in self.input_iedb_files:
            with open(input_iedb_file, 'r') as reader:
                iedb_tsv_reader = csv.DictReader(reader, delimiter='\t')
                filename = os.path.basename(input_iedb_file)

                pattern = re.compile(rf"{re.escape(self.sample_name)}\.(\w+(?:-\d+\.\d+)?)")
                match = pattern.match(filename)
                method = match.group(1)

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

                    if tsv_indices_from_label[protein_label] is not None:
                        tsv_indices = tsv_indices_from_label[protein_label]

                    for index in tsv_indices:
                        key = '|'.join([index, position])
                        if key not in iedb_results:
                            iedb_results[key]                      = {}
                            iedb_results[key]['mt_scores']         = {}
                            iedb_results[key]['mt_epitope_seq']    = epitope
                            iedb_results[key]['position']          = position
                            iedb_results[key]['tsv_index']         = index
                            iedb_results[key]['allele']            = allele
                        iedb_results[key]['mt_scores'].update(scores)
        return iedb_results

    def add_summary_metrics(self, iedb_results):
        iedb_results_with_metrics = {}
        for key, result in iedb_results.items():
            for metric in ['ic50', 'ic50_percentile', 'immunogenicity_percentile', 'presentation_percentile', 'percentile']:
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

                    if metric == 'ic50':
                        result['median_mt_{}'.format(metric)] = median([score['ic50'] for score in mt_values.values()])
                    else:
                        result['median_mt_{}'.format(metric)] = median([score['percentile'] for score in mt_values.values()])
                iedb_results_with_metrics[key]  = result
        return iedb_results_with_metrics

    def process_input_iedb_file(self):
        iedb_results              = self.parse_iedb_file()
        iedb_results_with_metrics = self.add_summary_metrics(iedb_results)
        return iedb_results_with_metrics

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
            'Median IC50 Percentile',
            'Best IC50 Percentile',
            'Best IC50 Percentile Method',
            'Median Immunogenicity Percentile',
            'Best Immunogenicity Percentile',
            'Best Immunogenicity Percentile Method',
            'Median Presentation Percentile',
            'Best Presentation Percentile',
            'Best Presentation Percentile Method',
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
            if method == 'MixMHCpred':
                headers.append("%s Binding Score" % pretty_method)
            elif method in ['BigMHC_EL', 'netmhciipan_el', 'netmhcpan_el']:
                headers.append("%s Presentation Score" % pretty_method)
            elif method in ['BigMHC_IM', 'DeepImmuno', 'PRIME']:
                headers.append("%s Immunogenicity Score" % pretty_method)
            else:
                headers.append("%s IC50 Score" % pretty_method)
            headers.append("%s Percentile" % pretty_method)
        if self.add_sample_name:
            headers.append("Sample Name")
        return headers

    def flurry_headers(self, headers):
        headers.append("MHCflurryEL Processing Score")
        headers.append("MHCflurryEL Processing Percentile")
        headers.append("MHCflurryEL Presentation Score")
        headers.append("MHCflurryEL Presentation Percentile")

    def add_prediction_scores(self, row, mt_scores):
        for method in self.prediction_methods():
            pretty_method = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)
            if pretty_method == 'MHCflurry':
                if self.flurry_state == 'EL_only' or self.flurry_state == 'both':
                    row['MHCflurryEL Processing Score'] = self.score_or_na(mt_scores, 'MHCflurryEL Processing', 'presentation')
                    row['MHCflurryEL Processing Percentile'] = self.score_or_na(mt_scores, 'MHCflurryEL Processing', 'percentile')
                    row['MHCflurryEL Presentation Score'] = self.score_or_na(mt_scores, 'MHCflurryEL Presentation', 'presentation')
                    row['MHCflurryEL Presentation Percentile'] = self.score_or_na(mt_scores, 'MHCflurryEL Presentation', 'percentile')
                if self.flurry_state in ['both', 'BA_only', None]:
                    row['MHCflurry IC50 Score'] = self.score_or_na(mt_scores, 'MHCflurry', 'ic50')
                    row['MHCflurry Percentile'] = self.score_or_na(mt_scores, 'MHCflurry', 'percentile')
            else:
                if pretty_method == 'MixMHCpred':
                    row[f'{pretty_method} Binding Score'] = self.score_or_na(mt_scores, pretty_method, 'binding_score')
                elif pretty_method in ['BigMHC_EL', 'NetMHCIIpanEL', 'NetMHCpanEL', 'MixMHCpred']:
                    row[f'{pretty_method} Presentation Score'] = self.score_or_na(mt_scores, pretty_method, 'presentation')
                elif pretty_method in ['BigMHC_IM', 'DeepImmuno', 'PRIME']:
                    row[f'{pretty_method} Immunogenicity Score'] = self.score_or_na(mt_scores, pretty_method, 'immunogenicity')
                else:
                    row[f'{pretty_method} IC50 Score'] = self.score_or_na(mt_scores, pretty_method, 'ic50')
                row[f'{pretty_method} Percentile'] = self.score_or_na(mt_scores, pretty_method, 'percentile')
        return row

    def execute(self):
        tmp_output_file = self.output_file + '.tmp'
        tmp_output_filehandle = open(tmp_output_file, 'w')
        tsv_writer = csv.DictWriter(tmp_output_filehandle, delimiter='\t', fieldnames=self.output_headers())
        tsv_writer.writeheader()

        iedb_results = self.process_input_iedb_file()
        for result in iedb_results.values():
            row = {
                'HLA Allele'          : result['allele'],
                'Sub-peptide Position': result['position'],
                'Epitope Seq'         : result['mt_epitope_seq'],
                'Mutation'            : result['tsv_index'],
                #Median IC50 Score
                'Median IC50 Score': self.rounded_score_or_na(result['median_mt_ic50']),
                #Median Percentile
                'Median Percentile': self.rounded_score_or_na(result['median_mt_percentile']),
                #Median IC50 Percentile
                'Median IC50 Percentile': self.rounded_score_or_na(result['median_mt_ic50_percentile']),
                #Median Immunogenicity Percentile
                'Median Immunogenicity Percentile': self.rounded_score_or_na(result['median_mt_immunogenicity_percentile']),
                #Median Presentation Percentile
                'Median Presentation Percentile': self.rounded_score_or_na(result['median_mt_presentation_percentile']),
                #Best IC50 Score
                'Best IC50 Score': self.rounded_score_or_na(result['best_mt_ic50']),
                'Best IC50 Score Method': result['best_mt_ic50_method'],
                #Best Percentile
                'Best Percentile': self.rounded_score_or_na(result['best_mt_percentile']),
                'Best Percentile Method': result['best_mt_percentile_method'],
                #Best IC50 Percentile
                'Best IC50 Percentile': self.rounded_score_or_na(result['best_mt_ic50_percentile']),
                'Best IC50 Percentile Method': result['best_mt_ic50_percentile_method'],
                #Best Immunogenicity Percentile
                'Best Immunogenicity Percentile': self.rounded_score_or_na(result['best_mt_immunogenicity_percentile']),
                'Best Immunogenicity Percentile Method': result['best_mt_immunogenicity_percentile_method'],
                #Best Presentation Percentile
                'Best Presentation Percentile': self.rounded_score_or_na(result['best_mt_presentation_percentile']),
                'Best Presentation Percentile Method': result['best_mt_presentation_percentile_method'],
            }
            row = self.add_prediction_scores(row, result['mt_scores'])
            if self.add_sample_name:
                row['Sample Name'] = self.sample_name
            tsv_writer.writerow(row)

        tmp_output_filehandle.close()
        os.replace(tmp_output_file, self.output_file)


class PvacspliceOutputParser(UnmatchedSequencesOutputParser):
    def parse_iedb_file(self):
        # input key file
        with open(self.key_file, 'r') as key_file_reader:
            protein_identifiers_from_label = yaml.load(key_file_reader, Loader=yaml.FullLoader)
        # final output
        iedb_results = {}
        for input_iedb_file in self.input_iedb_files:
            # input iedb file
            with open(input_iedb_file, 'r') as reader:
                iedb_tsv_reader = csv.DictReader(reader, delimiter='\t')
                filename = os.path.basename(input_iedb_file)
                pattern = re.compile(rf"{re.escape(self.sample_name)}\.(\w+(?:-\d+\.\d+)?)")
                match = pattern.match(filename)
                method = match.group(1)

                # header: allele, seq_num, start, end, length, peptide, ic50, percentile_rank
                for line in iedb_tsv_reader:
                    if "Warning: Potential DNA sequence(s)" in line['allele']:
                        continue
                    allele         = line['allele']
                    fasta_label    = int(line['seq_num'])
                    epitope        = line['peptide']
                    peptide_length = len(epitope)
                    scores         = self.get_scores(line, method)
                    # get fasta_id/combined_name from fasta key file
                    if protein_identifiers_from_label[fasta_label] is not None:
                        # comma-separated string (1 or more ids) as 1 entry in list
                        protein_label = protein_identifiers_from_label[fasta_label][0]
                        # one index at a time
                        for key in protein_label.split(','):

                            if key not in iedb_results:
                                iedb_results[key]                   = {}
                                iedb_results[key]['mt_scores']      = {}
                                iedb_results[key]['mt_epitope_seq'] = epitope
                                iedb_results[key]['fasta_id']       = fasta_label
                                iedb_results[key]['tsv_index']      = key
                                iedb_results[key]['allele']         = allele
                                iedb_results[key]['peptide_length'] = peptide_length
                            iedb_results[key]['mt_scores'].update(scores)

        return iedb_results

    def base_headers(self):
        return[
            'Chromosome',
            'Start',
            'Stop',
            'Reference',
            'Variant',
            'Junction',
            'Junction Start',
            'Junction Stop',
            'Junction Score',
            'Junction Anchor',
            'Transcript',
            'Transcript Support Level',
            'Canonical',
            'MANE Select',
            'Biotype',
            'Transcript CDS Flags',
            'Ensembl Gene ID',
            'Variant Type',
            'Amino Acid Change',
            'Gene Name',
            'HGVSc',
            'HGVSp',
            'WT Protein Length',
            'ALT Protein Length',
            'Frameshift Event',
            'Protein Position', # start position of peptide in alt protein
            'HLA Allele',
            'Peptide Length',
            'Epitope Seq',
            'Median IC50 Score',
            'Best IC50 Score',
            'Best IC50 Score Method',
            'Median Percentile',
            'Best Percentile',
            'Best Percentile Method',
            'Median IC50 Percentile',
            'Best IC50 Percentile',
            'Best IC50 Percentile Method',
            'Median Immunogenicity Percentile',
            'Best Immunogenicity Percentile',
            'Best Immunogenicity Percentile Method',
            'Median Presentation Percentile',
            'Best Presentation Percentile',
            'Best Presentation Percentile Method',
            'Tumor DNA Depth',
            'Tumor DNA VAF',
            'Tumor RNA Depth',
            'Tumor RNA VAF',
            'Normal Depth',
            'Normal VAF',
            'Gene Expression',
            'Transcript Expression',
            'Index', # this is junction index
            'Fasta Key', # unique num for traceback to correct sequence - key to combined fasta header
        ]

    def execute(self):
        tmp_output_file = self.output_file + '.tmp'
        tmp_output_filehandle = open(tmp_output_file, 'w')
        tsv_writer = csv.DictWriter(tmp_output_filehandle, delimiter='\t', fieldnames=self.output_headers())
        tsv_writer.writeheader()

        # added for pvacsplice - variant info
        tsv_entries = self.parse_input_tsv_file()

        # get binding info from iedb files
        iedb_results = self.process_input_iedb_file()

        # from input iedb files
        for result in iedb_results.values():
            # get unique index
            (final_index, protein_position) = result['tsv_index'].rsplit('.', 1)
            tsv_entry = tsv_entries[final_index]
            row = {
                'Chromosome'          : tsv_entry['chromosome_name'],
                'Start'               : tsv_entry['start'],
                'Stop'                : tsv_entry['stop'],
                'Reference'           : tsv_entry['reference'],
                'Variant'             : tsv_entry['variant'],
                'Transcript'          : tsv_entry['transcript_name'],
                'Transcript Support Level': tsv_entry['transcript_support_level'],
                'Canonical'           : tsv_entry['canonical'],
                'MANE Select'         : tsv_entry['mane_select'],
                'Biotype'             : tsv_entry['biotype'],
                'Transcript CDS Flags': tsv_entry['transcript_cds_flags'],
                ### junction info from RegTools
                'Junction'            : tsv_entry['name'],
                'Junction Start'      : tsv_entry['junction_start'],
                'Junction Stop'       : tsv_entry['junction_stop'],
                'Junction Score'      : tsv_entry['score'],
                'Junction Anchor'     : tsv_entry['anchor'],
                ###
                'Ensembl Gene ID'     : tsv_entry['gene_name'],
                'Variant Type'        : tsv_entry['variant_type'],
                'Amino Acid Change'   : tsv_entry['amino_acid_change'],
                'Protein Position' : protein_position,
                'Gene Name'           : tsv_entry['gene_name'],
                'HGVSc'               : tsv_entry['hgvsc'],
                'HGVSp'               : tsv_entry['hgvsp'],
                'Index'               : final_index,
                'Fasta Key'           : result['fasta_id'],
                'WT Protein Length' : tsv_entry['wt_protein_length'],
                'ALT Protein Length': tsv_entry['alt_protein_length'],
                'Frameshift Event'     : tsv_entry['frameshift_event'],
                ### pvacbind info
                'HLA Allele'          : result['allele'],
                'Peptide Length'      : len(result['mt_epitope_seq']),
                'Epitope Seq'         : result['mt_epitope_seq'],
                #Median IC50 Score
                'Median IC50 Score': self.rounded_score_or_na(result['median_mt_ic50']),
                #Median Percentile
                'Median Percentile': self.rounded_score_or_na(result['median_mt_percentile']),
                #Median IC50 Percentile
                'Median IC50 Percentile': self.rounded_score_or_na(result['median_mt_ic50_percentile']),
                #Median Immunogenicity Percentile
                'Median Immunogenicity Percentile': self.rounded_score_or_na(result['median_mt_immunogenicity_percentile']),
                #Median Presentation Percentile
                'Median Presentation Percentile': self.rounded_score_or_na(result['median_mt_presentation_percentile']),
                #Best IC50 Score
                'Best IC50 Score': self.rounded_score_or_na(result['best_mt_ic50']),
                'Best IC50 Score Method': result['best_mt_ic50_method'],
                #Best Percentile
                'Best Percentile': self.rounded_score_or_na(result['best_mt_percentile']),
                'Best Percentile Method': result['best_mt_percentile_method'],
                #Best IC50 Percentile
                'Best IC50 Percentile': self.rounded_score_or_na(result['best_mt_ic50_percentile']),
                'Best IC50 Percentile Method': result['best_mt_ic50_percentile_method'],
                #Best Immunogenicity Percentile
                'Best Immunogenicity Percentile': self.rounded_score_or_na(result['best_mt_immunogenicity_percentile']),
                'Best Immunogenicity Percentile Method': result['best_mt_immunogenicity_percentile_method'],
                #Best Presentation Percentile
                'Best Presentation Percentile': self.rounded_score_or_na(result['best_mt_presentation_percentile']),
                'Best Presentation Percentile Method': result['best_mt_presentation_percentile_method'],
            }
            row = self.add_prediction_scores(row, result['mt_scores'])

            for (tsv_key, row_key) in zip(['gene_expression', 'transcript_expression', 'normal_vaf', 'tdna_vaf', 'trna_vaf'], ['Gene Expression', 'Transcript Expression', 'Normal VAF', 'Tumor DNA VAF', 'Tumor RNA VAF']):
                if tsv_key in tsv_entry:
                    if tsv_entry[tsv_key] == 'NA':
                        row[row_key] = 'NA'
                    else:
                        # no --normal-sample-name parameter causes ValueError here bc tries to convert empty string to float
                        if 'normal' in tsv_key and tsv_entry[tsv_key] == '':
                            row[row_key] = 'NA'
                        else:
                            row[row_key] = round(float(tsv_entry[tsv_key]), 3)

            for (tsv_key, row_key) in zip(['normal_depth', 'tdna_depth', 'trna_depth'], ['Normal Depth', 'Tumor DNA Depth', 'Tumor RNA Depth']):
                if tsv_key in tsv_entry:
                    row[row_key] = tsv_entry[tsv_key]
                elif 'normal' in tsv_key and tsv_entry[tsv_key] == '':
                    row[row_key] = 'NA'

            if self.add_sample_name:
                row['Sample Name'] = self.sample_name
            tsv_writer.writerow(row)

        tmp_output_filehandle.close()
        os.replace(tmp_output_file, self.output_file)
