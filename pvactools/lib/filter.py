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
                for criteria in self.filter_criteria:
                    value = line[criteria['column']]
                    if value == 'inf':
                        value = sys.maxsize
                    if value == 'NA':
                        if criteria['exclude_nas']:
                            to_filter = True
                    else:
                        if not eval("{} {} {}".format(value, criteria['operator'], criteria['threshold'])):
                            to_filter = True
                if not to_filter:
                    writer.writerow(line)
