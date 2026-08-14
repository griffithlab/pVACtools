import pandas as pd
import csv
import operator
import sys

pd.options.mode.chained_assignment = None

FILTER_OPERATORS = {
    "<": operator.lt,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne,
    ">=": operator.ge,
    ">": operator.gt,
}

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

                def process_criterion(criterion):
                    value = line[criterion.column]
                    if value == 'NA' or criterion.skip_value == value:
                        return False
                    if value == 'inf':
                        numeric_value = sys.maxsize
                    else:
                        try:
                            numeric_value = float(value)
                        except (TypeError, ValueError):
                            raise ValueError(
                                "Invalid numeric value {!r} for filter column {!r}.".format(
                                    value,
                                    criterion.column,
                                )
                            )
                    return not criterion.comparator(numeric_value, criterion.threshold)

                if self.filter_strategy == "AND":
                    to_filter = any(process_criterion(criterion) for criterion in self.filter_criteria)
                else:
                    to_filter = all(process_criterion(criterion) for criterion in self.filter_criteria)

                if not to_filter:
                    writer.writerow(line)

class FilterCriterion:
    def __init__(self, column, operator, threshold, skip_value=None):
        if operator not in FILTER_OPERATORS:
            raise ValueError(
                "Unsupported filter operator {!r}. Supported operators are: {}.".format(
                    operator,
                    ", ".join(FILTER_OPERATORS),
                )
            )
        try:
            numeric_threshold = float(threshold)
        except (TypeError, ValueError):
            raise ValueError(
                "Invalid numeric threshold {!r} for filter column {!r}.".format(
                    threshold,
                    column,
                )
            )
        self.column = column
        self.operator = operator
        self.comparator = FILTER_OPERATORS[operator]
        self.threshold = numeric_threshold
        self.skip_value = skip_value
