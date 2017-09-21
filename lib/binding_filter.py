import argparse
import sys
import re
import csv

class BindingFilter:
    def __init__(self, input_file, output_file, binding_threshold, minimum_fold_change, top_score_metric):
        self.input_file = input_file
        self.output_file = output_file
        self.binding_threshold = binding_threshold
        self.minimum_fold_change = minimum_fold_change
        self.top_score_metric = top_score_metric

    def execute(self):
        with open(self.input_file, 'r') as csv_input_file:
            reader = csv.DictReader(csv_input_file, delimiter='\t')
            fieldnames = reader.fieldnames

            with open(self.output_file, 'w') as csv_output_file:
                writer = csv.DictWriter(
                    csv_output_file,
                    fieldnames,
                    delimiter = '\t',
                    lineterminator = '\n'
                )
                writer.writeheader()

                for entry in reader:
                    if self.top_score_metric == 'median':
                        score = float(entry['Median MT Score'])
                        fold_change = sys.maxsize if entry['Median Fold Change'] == 'NA' else float(entry['Median Fold Change'])
                    elif self.top_score_metric == 'lowest':
                        score = float(entry['Best MT Score'])
                        fold_change = sys.maxsize if entry['Corresponding Fold Change'] == 'NA' else float(entry['Corresponding Fold Change'])

                    if score > self.binding_threshold or fold_change < self.minimum_fold_change:
                        continue

                    writer.writerow(entry)

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser('%s binding_filter' % tool)
        parser.add_argument(
            'input_file',
            help="The final report .tsv file to filter"
        )
        parser.add_argument(
            'output_file',
            help="Output .tsv file containing list of filtered "
                 + "epitopes based on binding affinity"
        )
        parser.add_argument(
            '-b', '--binding-threshold', type=int,
            help="Report only epitopes where the mutant allele "
                 + "has ic50 binding scores below this value. Default: 500",
            default=500
        )
        parser.add_argument(
            '-c', '--minimum-fold-change', type=int,
            help="Minimum fold change between mutant binding "
                 + "score and wild-type score. The default is 0, which "
                 + "filters no results, but 1 is often a sensible "
                 + "option (requiring that binding is better to the MT than WT). "
                 + "Default: 0",
            default=0
        )
        parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            help="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. "
                 + "lowest: Best MT Score/Corresponding Fold Change - lowest MT ic50 binding score/corresponding fold change of all chosen prediction methods. "
                 + "median: Median MT Score/Median Fold Change - median MT ic50 binding score/fold change of all chosen prediction methods. "
                 + "Default: median",
            default='median',
        )
        return parser
