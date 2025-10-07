import argparse
import sys
import re
import csv
from pvactools.lib.filter import Filter, FilterCriterion
from pvactools.lib.allele_specific_binding_filter import AlleleSpecificBindingFilter
from pvactools.lib.run_utils import *

class BindingFilter:
    def __init__(
            self, input_file, output_file, binding_threshold=500, minimum_fold_change=None, top_score_metric='median', exclude_nas=False, allele_specific_binding_thresholds=False,
            binding_percentile_threshold=2.0, immunogenicity_percentile_threshold=2.0, presentation_percentile_threshold=2.0, percentile_threshold_strategy='conservative', file_type='pVACseq'
        ):
        self.input_file = input_file
        self.output_file = output_file
        self.binding_threshold = binding_threshold
        self.binding_percentile_threshold = binding_percentile_threshold
        self.immunogenicity_percentile_threshold = immunogenicity_percentile_threshold
        self.presentation_percentile_threshold = presentation_percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.minimum_fold_change = minimum_fold_change
        self.top_score_metric = top_score_metric
        self.exclude_nas = exclude_nas
        self.allele_specific_cutoffs = allele_specific_binding_thresholds
        self.file_type = file_type

    def execute(self):
        filter_criteria = []

        if self.allele_specific_cutoffs:
            AlleleSpecificBindingFilter(
                self.input_file, self.output_file,
                binding_threshold=self.binding_threshold,
                fold_change=self.minimum_fold_change,
                top_score_metric=self.top_score_metric,
                exclude_nas=self.exclude_nas,
                binding_precentile_threshold=self.binding_precentile_threshold,
                immunogenicity_percentile_threshold=self.immunogenicity_percentile_threshold,
                presentation_percentile_threshold=self.presentation_percentile_threshold,
                percentile_threshold_strategy=self.percentile_threshold_strategy,
                file_type=self.file_type
            ).execute()
        else:
            if self.file_type in ['pVACbind', 'pVACfuse', 'pVACsplice']:
                if self.top_score_metric == 'median':
                    ic50_column = 'Median IC50 Score'
                    binding_percentile_column = 'Median IC50 Percentile'
                    immunogenicity_percentile_column = 'Median Immunogenicity Percentile'
                    presentation_percentile_column = 'Median Presentation Percentile'
                elif self.top_score_metric == 'lowest':
                    ic50_column = 'Best IC50 Score'
                    binding_percentile_column = 'Best IC50 Percentile'
                    immunogenicity_percentile_column = 'Best Immunogenicity Percentile'
                    presentation_percentile_column = 'Best Presentation Percentile'
            else:
                if self.top_score_metric == 'median':
                    ic50_column = 'Median MT IC50 Score'
                    binding_percentile_column = 'Median MT IC50 Percentile'
                    immunogenicity_percentile_column = 'Median MT Immunogenicity Percentile'
                    presentation_percentile_column = 'Median MT Presentation Percentile'
                elif self.top_score_metric == 'lowest':
                    ic50_column = 'Best MT IC50 Score'
                    binding_percentile_column = 'Best MT IC50 Percentile'
                    immunogenicity_percentile_column = 'Best MT Immunogenicity Percentile'
                    presentation_percentile_column = 'Best MT Presentation Percentile'
            filter_criteria.append(FilterCriterion(ic50_column, '<=', self.binding_threshold, exclude_nas=self.exclude_nas))
            filter_criteria.append(FilterCriterion(binding_percentile_column, '<=', self.binding_percentile_threshold, exclude_nas=False))
            filter_criteria.append(FilterCriterion(immunogenicity_percentile_column, '<=', self.immunogenicity_percentile_threshold, exclude_nas=False))
            filter_criteria.append(FilterCriterion(presentation_percentile_column, '<=', self.presentation_percentile_threshold, exclude_nas=False))

            if self.minimum_fold_change is not None:
                if self.top_score_metric == 'median':
                    column = 'Median Fold Change'
                elif self.top_score_metric == 'lowest':
                    column = 'Corresponding Fold Change'
                filter_criteria.append(FilterCriterion(column, '>=', self.minimum_fold_change, exclude_nas=self.exclude_nas))
            Filter(self.input_file, self.output_file, filter_criteria, [], "AND" if self.percentile_threshold_strategy == 'conservative' else "OR").execute()

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s binding_filter' % tool,
            description="Filter variants processed by IEDB by binding score.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="The all_epitopes.tsv or filtered.tsv pVACseq report file to filter."
        )
        parser.add_argument(
            'output_file',
            help="Output .tsv file containing list of filtered "
                 + "epitopes based on binding affinity."
        )
        parser.add_argument(
            '-b', '--binding-threshold', type=int,
            help="Report only epitopes where the mutant allele "
                 + "has ic50 binding scores below this value.",
            default=500
        )
        parser.add_argument(
            '--binding-percentile-threshold', type=float_range(0.0,100.0),
            default=2.0,
            help="Report only epitopes where the mutant allele "
                 +"has a binding percentile rank below this value."
        )
        parser.add_argument(
            '--immunogenicity-percentile-threshold', type=float_range(0.0,100.0),
            default=2.0,
            help="Report only epitopes where the mutant allele "
                 +"has a immunogenicity percentile rank below this value."
        )
        parser.add_argument(
            '--presentation-percentile-threshold', type=float_range(0.0,100.0),
            default=2.0,
            help="Report only epitopes where the mutant allele "
                 +"has a presentation percentile rank below this value."
        )
        parser.add_argument(
            '--percentile-threshold-strategy',
            choices=['conservative', 'exploratory'],
            help="Specify the candidate inclusion strategy. The 'conservative' option requires a candidate to pass BOTH the binding threshold and percentile threshold (default)."
                 + " The 'exploratory' option requires a candidate to pass EITHER the binding threshold or the percentile threshold.",
            default="conservative",
        )
        if tool == 'pvacseq':
            parser.add_argument(
                '-c', '--minimum-fold-change', type=int,
                help="Minimum fold change between mutant binding "
                     + "score and wild-type score. The default is 0, which "
                     + "filters no results, but 1 is often a sensible "
                     + "option (requiring that binding is better to the MT than WT).",
                default=0
            )
        parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            help="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. "
                 + "lowest: Use the Best MT IC50 Score, Corresponding Fold Change, and Best MT Percentile "
                 + "(i.e. use the lowest MT ic50 binding score, orresponding fold change of all chosen prediction methods, and lowest MT percentile). "
                 + "median: Use the Median MT IC50 Score, Median Fold Change, and Median MT Percentile "
                 + "i.e. use the median MT ic50 binding score, fold change, and MT percentile of all chosen prediction methods).",
            default='median',
        )
        parser.add_argument(
            '--exclude-NAs',
            help="Exclude NA values from the filtered output.",
            default=False,
            action='store_true'
        )
        parser.add_argument(
            '-a', '--allele-specific-binding-thresholds',
            help="Use allele-specific binding thresholds. To print the allele-specific binding thresholds run `%s allele_specific_cutoffs`. " % tool
                 + "If an allele does not have a special threshold value, the `--binding-threshold` value will be used.",
            default=False,
            action='store_true',
        )
        return parser
