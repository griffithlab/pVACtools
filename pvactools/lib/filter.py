import pandas as pd
import csv
import sys

pd.options.mode.chained_assignment = None

class Filter:
    def __init__(self, input_file, output_file, filter_criteria, int_filter_columns=[]):
        self.input_file = input_file
        self.output_file = output_file
        self.filter_criteria = filter_criteria
        self.int_filter_columns = int_filter_columns

    def execute(self):
        with open(self.input_file, 'r') as read_fh, open(self.output_file, 'w') as write_fh:
            reader = csv.DictReader(read_fh, delimiter="\t")
            writer = csv.DictWriter(write_fh, delimiter='\t', fieldnames = reader.fieldnames, lineterminator="\n")
            writer.writeheader()
            for line in reader:
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
                        if not eval("{} {} {}".format(value, criterion.operator, criterion.threshold)):
                            to_filter = True
                if not to_filter:
                    writer.writerow(line)

class FilterCriterion:
    def __init__(self, column, operator, threshold, exclude_nas=False, skip_value=None):
        self.column = column
        self.operator = operator
        self.threshold = threshold
        self.exclude_nas = exclude_nas
        self.skip_value = skip_value
