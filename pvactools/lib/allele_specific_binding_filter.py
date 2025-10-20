import csv
import sys
from pvactools.lib.prediction_class import PredictionClass

class AlleleSpecificBindingFilter:
    def __init__(
            self, input_file, output_file, default_threshold=500, minimum_fold_change=None, top_score_metric='median',
            binding_percentile_threshold=2.0, immunogenicity_percentile_threshold=2.0, presentation_percentile_threshold=2.0, percentile_threshold_strategy="conservative", file_type='pVACseq'
        ):
        self.input_file = input_file
        self.output_file = output_file
        self.default_threshold = default_threshold
        self.minimum_fold_change = minimum_fold_change
        self.top_score_metric = top_score_metric
        self.binding_percentile_threshold = binding_percentile_threshold
        self.immunogenicity_percentile_threshold = immunogenicity_percentile_threshold
        self.presentation_percentile_threshold = presentation_percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.file_type = file_type

    def execute(self):
        with open(self.input_file, 'r') as input_fh:
            reader = csv.DictReader(input_fh, delimiter='\t')
            fieldnames = reader.fieldnames
            output_fh = open(self.output_file, 'w')
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

                threshold = PredictionClass.cutoff_for_allele(entry['HLA Allele'])
                threshold = self.default_threshold if threshold is None else float(threshold)

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
            output_fh.close()
