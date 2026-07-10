from abc import ABCMeta, abstractmethod
import sys
import csv
import operator
import os
import pandas as pd
from math import ceil, inf
from statistics import median
import yaml

from pvactools.lib.normalized_percentile_calculator import NormalizedPercentileCalculator

csv.field_size_limit(sys.maxsize)

class OutputParser(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        self.prediction_files        = kwargs['prediction_files']
        self.tsv_file                = kwargs['tsv_file']
        self.key_files               = kwargs['key_files']
        self.output_file             = kwargs['output_file']
        self.sample_name             = kwargs['sample_name']
        self.add_sample_name         = kwargs.get('add_sample_name_column')
        self.flurry_state            = kwargs.get('flurry_state')
        self.use_normalized_percentiles = kwargs.get('use_normalized_percentiles', False)
        reference_scores_path        = kwargs.get('reference_scores_path', '/tmp')
        self.normalized_percentile_calculator = NormalizedPercentileCalculator(reference_scores_path = reference_scores_path)

    def parse_tsv_file(self):
        with open(self.tsv_file, 'r') as reader:
            tsv_reader = csv.DictReader(reader, delimiter='\t')
            tsv_entries = {}
            for line in tsv_reader:
                if line['index'] in tsv_entries:
                    sys.exit('Duplicate TSV indexes')
                tsv_entries[line['index']] = line
            return tsv_entries

    def process_prediction_files(self):
        prediction_results = self.parse_prediction_files()
        prediction_results_with_metrics = self.add_summary_metrics(prediction_results)
        return prediction_results_with_metrics

    def find_mutation_positions(self, wt_epitope_seq, mt_epitope_seq):
        mutated_positions = []
        for i,(wt_aa,mt_aa) in enumerate(zip(wt_epitope_seq,mt_epitope_seq)):
            if wt_aa != mt_aa:
                mutated_positions.append(i+1)
        if len(mutated_positions) == 0:
            return "NA"
        else:
            return ", ".join([str(x) for x in mutated_positions])

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
                percentile = self.normalized_percentile_calculator.calculate_normalized_percentile(
                    line.get('allele'),
                    len(line.get('peptide') or ''),
                    normalized_input,
                    method,
                    is_reversed,
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

        if m == 'immuscope_im':
            return self._make_score_entry(
                line, 'ImmuScope_IM', 'immunogenicity',
                line.get('ImmuScope_IM'), method,
                percentile_keys=None, percentile_fallback='NA', is_reversed=True
            )

        if m == 'netmhcpanel':
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
                entry['NetMHCpanEL']['percentile'] = self.normalized_percentile_calculator.calculate_normalized_percentile(
                    line.get('allele'),
                    len(line.get('peptide') or ''),
                    self._parse_float_or_na(presentation),
                    'NetMHCpanEL',
                    is_reversed=True
                )

            return entry

        if 'netmhciipanel' in m:
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
                entry['NetMHCIIpanEL']['percentile'] = self.normalized_percentile_calculator.calculate_normalized_percentile(
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

        if m == 'mixmhc2pred':
            return self._make_score_entry(
                line, 'MixMHC2pred', 'presentation',
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

        percentile_keys = [
            key for key in ('percentile', 'percentile_rank', 'rank') if key in line
        ] or None

        return self._make_score_entry(
            line, method, 'ic50',
            line.get('ic50'), method,
            percentile_keys=percentile_keys
        )

    def format_match_na(self, result, metric):
        return {method: {field: 'NA' for field in fields.keys()} for method, fields in result[f'mt_{metric}s'].items()}

    def get_values_for_summary_metrics(self, result, metric, epitope_type):
        metric_values = dict()
        if metric in ['ic50', 'percentile']:
            for (method, values) in result['{}_scores'.format(epitope_type)].items():
                metric_values[method] = {field: score for field, score in values.items() if field == metric and score != 'NA'}
                if not metric_values[method]:
                    del metric_values[method]
        else:
            if metric == 'ic50_percentile':
                #also include MixMHCpred percentiles in ic50_percentile
                result_subset = {method: scores for method, scores in result['{}_scores'.format(epitope_type)].items() if 'ic50' in scores or 'binding_score' in scores}
            else:
                field = metric.replace("_percentile", "")
                result_subset = {method: scores for method, scores in result['{}_scores'.format(epitope_type)].items() if field in scores}
            for (method, values) in result_subset.items():
                metric_values[method] = {field: score for field, score in values.items() if field == 'percentile' and score != 'NA'}
                if not metric_values[method]:
                    del metric_values[method]
        return metric_values

    def prediction_methods(self):
        methods = set()
        for prediction_file in self.prediction_files:
            filename = os.path.basename(prediction_file)
            method = filename.rsplit('.', 5)[1]
            methods.add(method)

        return sorted(list(methods))

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


class MatchedSequencesOutputParser(OutputParser):
    def execute(self):
        tmp_output_file = self.output_file + '.tmp'
        tmp_output_filehandle = open(tmp_output_file, 'w')
        tsv_writer = csv.DictWriter(tmp_output_filehandle, delimiter='\t', fieldnames=self.output_headers())
        tsv_writer.writeheader()

        tsv_entries = self.parse_tsv_file()

        prediction_results = self.process_prediction_files()

        for result in prediction_results.values():
            tsv_index = result['tsv_index']
            tsv_entry = tsv_entries[tsv_index]

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

            row = self.construct_row(result, tsv_entry, corresponding_fold_change, median_fold_change)
            tsv_writer.writerow(row)

        tmp_output_filehandle.close()
        os.replace(tmp_output_file, self.output_file)


    def parse_prediction_files(self):
        # input key file
        protein_identifiers_from_label = {}
        for key_file in self.key_files:
            with open(key_file, 'r') as key_file_reader:
                chunk = key_file.rsplit('.', 2)[1]
                protein_identifiers_from_label[chunk] = yaml.load(key_file_reader, Loader=yaml.FullLoader)
        # final output
        prediction_results = {}
        wt_prediction_results = {}
        for prediction_file in self.prediction_files:
            with open(prediction_file, 'r') as reader:
                chunk = prediction_file.rsplit('.', 2)[1]
                prediction_reader = csv.DictReader(reader, delimiter='\t')
                filename = os.path.basename(prediction_file)
                method = filename.rsplit('.', 5)[1]

                # header: allele, seq_num, start, end, length, peptide, ic50, percentile_rank
                for line in prediction_reader:
                    if "Warning: Potential DNA sequence(s)" in line['allele']:
                        continue
                    allele         = line['allele']
                    fasta_label    = int(line['seq_num'])
                    epitope        = line['peptide']
                    peptide_length = len(epitope)
                    scores         = self.get_scores(line, method)
                    # get fasta_id/combined_name from fasta key file
                    if protein_identifiers_from_label[chunk][fasta_label] is not None:
                        # comma-separated string (1 or more ids) as 1 entry in list
                        protein_labels = protein_identifiers_from_label[chunk][fasta_label]
                        # one index at a time
                        for key in protein_labels:
                            (protein_type, rest) = key.split('.', 1)
                            (tsv_index, position) = rest.rsplit('|', 1)
                            if protein_type in ['ALT', 'MT']:
                                if rest not in prediction_results:
                                    prediction_results[rest]                   = {}
                                    prediction_results[rest]['mt_scores']      = {}
                                    prediction_results[rest]['mt_epitope_seq'] = epitope
                                    prediction_results[rest]['fasta_id']       = fasta_label
                                    prediction_results[rest]['tsv_index']      = tsv_index
                                    prediction_results[rest]['allele']         = allele
                                    prediction_results[rest]['peptide_length'] = peptide_length
                                    prediction_results[rest]['position']       = int(position) + 1
                                prediction_results[rest]['mt_scores'].update(scores)
                            else:
                                if tsv_index not in wt_prediction_results:
                                    wt_prediction_results[tsv_index] = {}
                                if position not in wt_prediction_results[tsv_index]:
                                    wt_prediction_results[tsv_index][position] = {}
                                    wt_prediction_results[tsv_index][position]['wt_scores'] = {}
                                wt_prediction_results[tsv_index][position]['wt_epitope_seq'] = epitope
                                wt_prediction_results[tsv_index][position]['wt_scores'].update(scores)

        return self.match_wildtype_and_mutant_entries(prediction_results, wt_prediction_results)

    def match_wildtype_and_mutant_entries(self, prediction_results, wt_prediction_results):
        for key, mt_result in prediction_results.items():
            (tsv_index, position) = key.rsplit('|', 1)
            if tsv_index in wt_prediction_results and position in wt_prediction_results[tsv_index]:
                wt_result = wt_prediction_results[tsv_index][position]
                mt_result['wt_epitope_seq'] = wt_result['wt_epitope_seq']
                mt_result['wt_scores']      = wt_result['wt_scores']
                mt_result['mutation_position'] = self.find_mutation_positions(wt_result['wt_epitope_seq'], mt_result['mt_epitope_seq'])
            else:
                mt_result['wt_epitope_seq'] = 'NA'
                mt_result['wt_scores']      = self.format_match_na(mt_result, 'score')
                mt_result['mutation_position'] = 'NA'
        return prediction_results

    def add_prediction_scores(self, row, mt_scores, wt_scores):
        for method in self.prediction_methods():
            if method == 'MHCflurry':
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
                if method in ['MixMHCpred']:
                    row[f'{method} MT Binding Score'] = self.score_or_na(mt_scores, method, 'binding_score')
                    row[f'{method} WT Binding Score'] = self.score_or_na(wt_scores, method, 'binding_score')
                elif method in ['BigMHC_EL', 'NetMHCIIpanEL', 'NetMHCpanEL', 'MixMHC2pred']:
                    row[f'{method} MT Presentation Score'] = self.score_or_na(mt_scores, method, 'presentation')
                    row[f'{method} WT Presentation Score'] = self.score_or_na(wt_scores, method, 'presentation')
                elif method in ['BigMHC_IM', 'DeepImmuno', 'PRIME', 'ImmuScope_IM']:
                    row[f'{method} MT Immunogenicity Score'] = self.score_or_na(mt_scores, method, 'immunogenicity')
                    row[f'{method} WT Immunogenicity Score'] = self.score_or_na(wt_scores, method, 'immunogenicity')
                else:
                    row[f'{method} MT IC50 Score'] = self.score_or_na(mt_scores, method, 'ic50')
                    row[f'{method} WT IC50 Score'] = self.score_or_na(wt_scores, method, 'ic50')
                row[f'{method} MT Percentile'] = self.score_or_na(mt_scores, method, 'percentile')
                row[f'{method} WT Percentile'] = self.score_or_na(wt_scores, method, 'percentile')
        return row

    def output_headers(self):
        headers = self.base_headers()
        for method in self.prediction_methods():
            if method.lower() == 'mhcflurry':
                if self.flurry_state == 'EL_only':
                    self.flurry_headers(headers)
                    continue
                elif self.flurry_state == 'both':
                    self.flurry_headers(headers)

            if method in ['MixMHCpred']:
                headers.append("%s WT Binding Score" % method)
                headers.append("%s MT Binding Score" % method)
            elif method in ['BigMHC_EL', 'NetMHCIIpanEL', 'NetMHCpanEL', 'MixMHC2pred']:
                headers.append("%s WT Presentation Score" % method)
                headers.append("%s MT Presentation Score" % method)
            elif method in ['BigMHC_IM', 'DeepImmuno', 'PRIME', 'ImmuScope_IM']:
                headers.append("%s WT Immunogenicity Score" % method)
                headers.append("%s MT Immunogenicity Score" % method)
            else:
                headers.append("%s WT IC50 Score" % method)
                headers.append("%s MT IC50 Score" % method)
            headers.append("%s WT Percentile" % method)
            headers.append("%s MT Percentile" % method)
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

    def add_summary_metrics(self, prediction_results):
        prediction_results_with_metrics = {}
        for key, result in prediction_results.items():
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

                prediction_results_with_metrics[key]  = result

        return prediction_results_with_metrics

class UnmatchedSequencesOutputParser(OutputParser):
    def execute(self):
        tmp_output_file = self.output_file + '.tmp'
        tmp_output_filehandle = open(tmp_output_file, 'w')
        tsv_writer = csv.DictWriter(tmp_output_filehandle, delimiter='\t', fieldnames=self.output_headers())
        tsv_writer.writeheader()

        prediction_results = self.process_prediction_files()

        for result in prediction_results.values():
            row = self.construct_row(result)
            tsv_writer.writerow(row)

        tmp_output_filehandle.close()
        os.replace(tmp_output_file, self.output_file)

    def parse_prediction_files(self):
        protein_identifiers_from_label = {}
        for key_file in self.key_files:
            with open(key_file, 'r') as key_file_reader:
                chunk = key_file.rsplit('.', 2)[1]
                protein_identifiers_from_label[chunk] = yaml.load(key_file_reader, Loader=yaml.FullLoader)
        prediction_results = {}
        for prediction_file in self.prediction_files:
            with open(prediction_file, 'r') as reader:
                chunk = prediction_file.rsplit('.', 2)[1]
                prediction_reader = csv.DictReader(reader, delimiter='\t')
                filename = os.path.basename(prediction_file)
                method = filename.rsplit('.', 5)[1]

                for line in prediction_reader:
                    if "Warning: Potential DNA sequence(s)" in line['allele']:
                        continue
                    fasta_label  = int(line['seq_num'])
                    epitope        = line['peptide']
                    scores         = self.get_scores(line, method)
                    allele         = line['allele']
                    peptide_length = len(epitope)

                    if protein_identifiers_from_label[chunk][fasta_label] is not None:
                        protein_labels = protein_identifiers_from_label[chunk][fasta_label]

                    for key in protein_labels:
                        (tsv_index, position) = key.rsplit('|', 1)
                        if 'core_peptide' in line and int(line['end']) - int(line['start']) == 8:
                            #Start and end refer to the position of the core peptide
                            #Infer the (start) position of the peptide from the positions of the core peptide
                            position = int(position) - line['peptide'].find(line['core_peptide'])

                        if key not in prediction_results:
                            prediction_results[key]                      = {}
                            prediction_results[key]['mt_scores']         = {}
                            prediction_results[key]['mt_epitope_seq']    = epitope
                            prediction_results[key]['position']          = int(position) + 1
                            prediction_results[key]['tsv_index']         = tsv_index
                            prediction_results[key]['allele']            = allele
                        prediction_results[key]['mt_scores'].update(scores)
        return prediction_results

    def add_summary_metrics(self, prediction_results):
        prediction_results_with_metrics = {}
        for key, result in prediction_results.items():
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
                prediction_results_with_metrics[key]  = result
        return prediction_results_with_metrics

    def output_headers(self):
        headers = self.base_headers()
        for method in self.prediction_methods():
            if method.lower() == 'mhcflurry':
                if self.flurry_state == 'EL_only':
                    self.flurry_headers(headers)
                    continue
                elif self.flurry_state == 'both':
                    self.flurry_headers(headers)

            if method in ['MixMHCpred']:
                headers.append("%s Binding Score" % method)
            elif method in ['BigMHC_EL', 'NetMHCIIpanEL', 'NetMHCpanEL', 'MixMHC2pred']:
                headers.append("%s Presentation Score" % method)
            elif method in ['BigMHC_IM', 'DeepImmuno', 'PRIME', 'ImmuScope_IM']:
                headers.append("%s Immunogenicity Score" % method)
            else:
                headers.append("%s IC50 Score" % method)
            headers.append("%s Percentile" % method)
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
            if method == 'MHCflurry':
                if self.flurry_state == 'EL_only' or self.flurry_state == 'both':
                    row['MHCflurryEL Processing Score'] = self.score_or_na(mt_scores, 'MHCflurryEL Processing', 'presentation')
                    row['MHCflurryEL Processing Percentile'] = self.score_or_na(mt_scores, 'MHCflurryEL Processing', 'percentile')
                    row['MHCflurryEL Presentation Score'] = self.score_or_na(mt_scores, 'MHCflurryEL Presentation', 'presentation')
                    row['MHCflurryEL Presentation Percentile'] = self.score_or_na(mt_scores, 'MHCflurryEL Presentation', 'percentile')
                if self.flurry_state in ['both', 'BA_only', None]:
                    row['MHCflurry IC50 Score'] = self.score_or_na(mt_scores, 'MHCflurry', 'ic50')
                    row['MHCflurry Percentile'] = self.score_or_na(mt_scores, 'MHCflurry', 'percentile')
            else:
                if method in ['MixMHCpred']:
                    row[f'{method} Binding Score'] = self.score_or_na(mt_scores, method, 'binding_score')
                elif method in ['BigMHC_EL', 'NetMHCIIpanEL', 'NetMHCpanEL', 'MixMHC2pred']:
                    row[f'{method} Presentation Score'] = self.score_or_na(mt_scores, method, 'presentation')
                elif method in ['BigMHC_IM', 'DeepImmuno', 'PRIME', 'ImmuScope_IM']:
                    row[f'{method} Immunogenicity Score'] = self.score_or_na(mt_scores, method, 'immunogenicity')
                else:
                    row[f'{method} IC50 Score'] = self.score_or_na(mt_scores, method, 'ic50')
                row[f'{method} Percentile'] = self.score_or_na(mt_scores, method, 'percentile')
        return row


class PvacbindOutputParser(UnmatchedSequencesOutputParser):
    def base_headers(self):
        return[
            'Index',
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

    def construct_row(self, result):
        row = {
            'HLA Allele'          : result['allele'],
            'Sub-peptide Position': result['position'],
            'Epitope Seq'         : result['mt_epitope_seq'],
            'Index'               : result['tsv_index'],
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
        return row

class PvacspliceOutputParser(MatchedSequencesOutputParser):
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
            'Junction Type',
            'Transcript',
            'Transcript Support Level',
            'Transcript Length',
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
            'Mutation Position',
            'HLA Allele',
            'Peptide Length',
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
            'Fasta Key', # unique num for traceback to correct sequence - key to combined fasta header
        ]

    def construct_row(self, result, tsv_entry, corresponding_fold_change, median_fold_change):
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
            'Transcript Length'   : tsv_entry['transcript_length'],
            ### junction info from RegTools
            'Junction'            : tsv_entry['name'],
            'Junction Start'      : tsv_entry['junction_start'],
            'Junction Stop'       : tsv_entry['junction_stop'],
            'Junction Score'      : tsv_entry['score'],
            'Junction Anchor'     : tsv_entry['anchor'],
            'Junction Type'       : tsv_entry['junction_type'],
            ###
            'Ensembl Gene ID'     : tsv_entry['gene_name'],
            'Variant Type'        : tsv_entry['variant_type'],
            'Amino Acid Change'   : tsv_entry['amino_acid_change'],
            'Protein Position'    : result['position'],
            'Gene Name'           : tsv_entry['gene_name'],
            'HGVSc'               : tsv_entry['hgvsc'],
            'HGVSp'               : tsv_entry['hgvsp'],
            'Index'               : result['tsv_index'],
            'Fasta Key'           : result['fasta_id'],
            'WT Protein Length'   : tsv_entry['wt_protein_length'],
            'ALT Protein Length'  : tsv_entry['alt_protein_length'],
            'Frameshift Event'    : tsv_entry['frameshift_event'],
            ### pvacbind info
            'HLA Allele'          : result['allele'],
            'Peptide Length'      : len(result['mt_epitope_seq']),
            'MT Epitope Seq'      : result['mt_epitope_seq'],
            'WT Epitope Seq'      : result['wt_epitope_seq'],
            'Mutation Position'   : result['mutation_position'] if 'mutation_position' in result else 'NA',
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
        return row

class PvacfuseOutputParser(MatchedSequencesOutputParser):
    def base_headers(self):
        return[
            'Chromosome',
            'Start',
            'Stop',
            'Transcript',
            'Gene Name',
            'Variant Type',
            'Read Support',
            'Expression',
            'Index',
            'HLA Allele',
            'Sub-peptide Position',
            'MT Epitope Seq',
            'WT Epitope Seq',
            'Mutation Position',
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

    def construct_row(self, result, tsv_entry, corresponding_fold_change, median_fold_change):
        row = {
            'Chromosome'          : tsv_entry['chromosome_name'],
            'Start'               : tsv_entry['start'],
            'Stop'                : tsv_entry['stop'],
            'Transcript'          : tsv_entry['transcript_name'],
            'Gene Name'           : tsv_entry['gene_name'],
            'Variant Type'        : tsv_entry['variant_type'],
            'Read Support'        : tsv_entry['fusion_read_support'],
            'Expression'          : tsv_entry['fusion_expression'],
            'Index'               : result['tsv_index'],
            'HLA Allele'          : result['allele'],
            'MT Epitope Seq'      : result['mt_epitope_seq'],
            'WT Epitope Seq'      : result['wt_epitope_seq'],
            'Mutation Position'   : result['mutation_position'] if 'mutation_position' in result else 'NA',
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

        if self.add_sample_name:
            row['Sample Name'] = self.sample_name

        return row

class PvacseqOutputParser(MatchedSequencesOutputParser):
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

    def construct_row(self, result, tsv_entry, corresponding_fold_change, median_fold_change):
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
            'Mutation'            : tsv_entry['amino_acid_change'],
            'Protein Position'    : tsv_entry['protein_position'],
            'Gene Name'           : tsv_entry['gene_name'],
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

        return row
