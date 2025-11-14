import argparse
import sys
import re
import csv

from pvactools.lib.prediction_class import PredictionClass
from pvactools.lib.filter import Filter, FilterCriterion
from pvactools.lib.run_utils import *

class BindingFilter:
    def __init__(
            self, input_file, output_file, binding_threshold=500, minimum_fold_change=None, top_score_metric='median', allele_specific_binding_thresholds=False,
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
        self.allele_specific_cutoffs = allele_specific_binding_thresholds
        self.hla_types = pd.read_csv(self.input_file, delimiter="\t", usecols=["HLA Allele"])['HLA Allele'].unique()
        per_allele_binding_thresholds = {}
        for hla_type in self.hla_types:
            threshold = PredictionClass.cutoff_for_allele(hla_type)
            if self.allele_specific_cutoffs and threshold is not None:
                per_allele_binding_thresholds[hla_type] = float(threshold)
            else:
                per_allele_binding_thresholds[hla_type] = binding_threshold
        self.per_allele_binding_thresholds = per_allele_binding_thresholds
        self.file_type = file_type

    def execute(self):
        with open(self.input_file, 'r') as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter='\t')
            fieldnames = reader.fieldnames
            writer = csv.DictWriter(output_fh, fieldnames, delimiter = '\t', lineterminator = '\n')
            writer.writeheader()

            for entry in reader:
                if self.file_type in ['pVACbind', 'pVACfuse', 'pVACsplice']:
                    if self.top_score_metric == 'median':
                        score = entry['Median IC50 Score']
                        binding_percentile = entry['Median IC50 Percentile']
                        immunogenicity_percentile = entry['Median Immunogenicity Percentile']
                        presentation_percentile = entry['Median Presentation Percentile']
                    elif self.top_score_metric == 'lowest':
                        score = entry['Best IC50 Score']
                        binding_percentile = entry['Best IC50 Percentile']
                        immunogenicity_percentile = entry['Best Immunogenicity Percentile']
                        presentation_percentile = entry['Best Presentation Percentile']
                else:
                    if self.top_score_metric == 'median':
                        score = entry['Median MT IC50 Score']
                        fold_change = sys.maxsize if entry['Median Fold Change'] == 'NA' else float(entry['Median Fold Change'])
                        binding_percentile = entry['Median MT IC50 Percentile']
                        immunogenicity_percentile = entry['Median MT Immunogenicity Percentile']
                        presentation_percentile = entry['Median MT Presentation Percentile']
                    elif self.top_score_metric == 'lowest':
                        score = entry['Best MT IC50 Score']
                        fold_change = sys.maxsize if entry['Corresponding Fold Change'] == 'NA' else float(entry['Corresponding Fold Change'])
                        binding_percentile = entry['Best MT IC50 Percentile']
                        immunogenicity_percentile = entry['Best MT Immunogenicity Percentile']
                        presentation_percentile = entry['Best MT Presentation Percentile']

                threshold = self.per_allele_binding_thresholds[entry['HLA Allele']]

                filters = [
                    False if score == 'NA' else float(score) > threshold,
                    False if binding_percentile == 'NA' else float(binding_percentile) > self.binding_percentile_threshold,
                    False if immunogenicity_percentile == 'NA' else float(immunogenicity_percentile) > self.immunogenicity_percentile_threshold,
                    False if presentation_percentile == 'NA' else float(presentation_percentile) > self.presentation_percentile_threshold,
                ]

                if self.percentile_threshold_strategy == 'conservative':
                    if any(filters):
                        continue
                elif self.percentile_threshold_strategy == 'exploratory':
                    if all(filters):
                        continue

                if self.minimum_fold_change is not None and fold_change < self.minimum_fold_change:
                    continue

                writer.writerow(entry)

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
            '-a', '--allele-specific-binding-thresholds',
            help="Use allele-specific binding thresholds. To print the allele-specific binding thresholds run `%s allele_specific_cutoffs`. " % tool
                 + "If an allele does not have a special threshold value, the `--binding-threshold` value will be used.",
            default=False,
            action='store_true',
        )
        return parser
