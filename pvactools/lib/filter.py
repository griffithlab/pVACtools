import pandas as pd
import csv
import sys

pd.options.mode.chained_assignment = None

class Filter:
    def __init__(self, input_file, output_file, filter_criteria, int_filter_columns=[], filter_strategy="AND"):
        self.input_file = input_file
        self.output_file = output_file
        self.filter_criteria = filter_criteria
        self.int_filter_columns = int_filter_columns
        self.filter_strategy = filter_strategy

    def execute(self):
        with open(self.input_file, 'r') as read_fh, open(self.output_file, 'w') as write_fh:
            reader = csv.DictReader(read_fh, delimiter="\t")
            writer = csv.DictWriter(write_fh, delimiter='\t', fieldnames = reader.fieldnames, lineterminator="\n")
            writer.writeheader()
            for line in reader:
                to_filter = False

                process_criterion = lambda criterion: (
                    (line[criterion.column] == 'NA' and criterion.exclude_nas) or
                    (line[criterion.column] != 'NA' and criterion.skip_value != line[criterion.column] and not eval("{} {} {}".format(
                        line[criterion.column] if line[criterion.column] != 'inf' else sys.maxsize,
                        criterion.operator,
                        criterion.threshold
                    )))
                )

                if self.filter_strategy == "AND":
                    to_filter = any(process_criterion(criterion) for criterion in self.filter_criteria)
                else:
                    to_filter = all(process_criterion(criterion) for criterion in self.filter_criteria)
                
                if not to_filter:
                    writer.writerow(line)

class FilterCriterion:
    def __init__(self, column, operator, threshold, exclude_nas=False, skip_value=None):
        self.column = column
        self.operator = operator
        self.threshold = threshold
        self.exclude_nas = exclude_nas
        self.skip_value = skip_value
