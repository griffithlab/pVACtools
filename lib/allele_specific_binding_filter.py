import csv
import sys
from lib.prediction_class import *

class AlleleSpecificBindingFilter:
    def __init__(self, input_file, output_file, default_threshold, minimum_fold_change, top_score_metric, exclude_nas):
        self.input_file = input_file
        self.output_file = output_file
        self.default_threshold = default_threshold
        self.minimum_fold_change = minimum_fold_change
        self.top_score_metric = top_score_metric
        self.exclude_nas = exclude_nas

    def execute(self):
        with open(self.input_file, 'r') as input_fh:
            reader = csv.DictReader(input_fh, delimiter='\t')
            fieldnames = reader.fieldnames
            output_fh = open(self.output_file, 'w')
            writer = csv.DictWriter(output_fh, fieldnames, delimiter = '\t', lineterminator = '\n')
            writer.writeheader()

            for entry in reader:
                if self.top_score_metric == 'median':
                    score = float(entry['Median MT Score'])
                    if self.exclude_nas and entry['Median Fold Change'] == 'NA':
                        continue
                    else:
                        fold_change = sys.maxsize if entry['Median Fold Change'] == 'NA' else float(entry['Median Fold Change'])
                elif self.top_score_metric == 'lowest':
                    score = float(entry['Best MT Score'])
                    if self.exclude_nas and entry['Corresponding Fold Change'] == 'NA':
                        continue
                    else:
                        fold_change = sys.maxsize if entry['Corresponding Fold Change'] == 'NA' else float(entry['Corresponding Fold Change'])

                threshold = PredictionClass.cutoff_for_allele(entry['HLA Allele'])
                threshold = self.default_threshold if threshold is None else float(threshold)

                if score > threshold or fold_change < self.minimum_fold_change:
                    continue

                writer.writerow(entry)
            output_fh.close()
