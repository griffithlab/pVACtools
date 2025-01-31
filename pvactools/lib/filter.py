import pandas as pd
import csv
import sys

pd.options.mode.chained_assignment = None

class Filter:
    def __init__(self, input_file, output_file, filter_criteria, int_filter_columns=[], percentile_threshold_strategy="conservative"):
        self.input_file = input_file
        self.output_file = output_file
        self.filter_criteria = filter_criteria
        self.int_filter_columns = int_filter_columns
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.percentile_columns = ['Median Percentile', 'Best Percentile', 'Median MT Percentile', 'Best MT Percentile']
        self.ic50_columns = ['Median IC50 Score', 'Best IC50 Score', 'Median MT IC50 Score', 'Best MT IC50 Score']

    def execute(self):
        with open(self.input_file, 'r') as read_fh, open(self.output_file, 'w') as write_fh:
            reader = csv.DictReader(read_fh, delimiter="\t")
            writer = csv.DictWriter(write_fh, delimiter='\t', fieldnames = reader.fieldnames, lineterminator="\n")
            writer.writeheader()
            for line in reader:
                binding_threshold_passed = False
                percentile_threshold_passed = False
                percentile_threshold_present = False
                to_filter = False

                for criterion in self.filter_criteria:
                    value = line[criterion.column]
                    if value == 'inf':
                        value = sys.maxsize
                    if value == 'NA':
                        if criterion.exclude_nas:
                            to_filter = True
                    elif criterion.skip_value is not None and value == criterion.skip_value:
                        pass
                    else:
                        condition_passes = eval("{} {} {}".format(value, criterion.operator, criterion.threshold))
                        if criterion.column in self.ic50_columns:
                            binding_threshold_passed = condition_passes
                        elif criterion.column in self.percentile_columns:
                            percentile_threshold_present = True
                            percentile_threshold_passed = condition_passes
                        elif not condition_passes:
                            to_filter = True

                if not to_filter:
                    if percentile_threshold_present:
                        if self.percentile_threshold_strategy == "conservative":
                            to_filter = not (binding_threshold_passed and percentile_threshold_passed)
                        else:
                            to_filter = not (binding_threshold_passed or percentile_threshold_passed)
                    else:
                        to_filter = not binding_threshold_passed

                if not to_filter:
                    writer.writerow(line)

class FilterCriterion:
    def __init__(self, column, operator, threshold, exclude_nas=False, skip_value=None):
        self.column = column
        self.operator = operator
        self.threshold = threshold
        self.exclude_nas = exclude_nas
        self.skip_value = skip_value
