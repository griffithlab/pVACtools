import pandas as pd

pd.options.mode.chained_assignment = None

class Filter:
    def __init__(self, input_file, output_file, filter_criteria, exclude_nas, int_filter_columns=[]):
        self.input_file = input_file
        self.output_file = output_file
        self.filter_criteria = filter_criteria
        self.exclude_nas = exclude_nas
        self.int_filter_columns = int_filter_columns

    def execute(self):
        data = pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False)
        header = data.columns
        clean_header = header.map(lambda x: x.replace(' ', '_') if isinstance(x, str) else x)
        data.columns = clean_header
        for criteria in self.filter_criteria:
            clean_column = criteria['column'].replace(' ', '_')
#           #clean_column != clean_column is a hacky way to keep all NA values
            if self.exclude_nas:
                expression = "(%s %s %s)" % (clean_column, criteria['operator'], criteria['threshold'])
            else:
                expression = "(%s %s %s) | (%s != %s)" % (clean_column, criteria['operator'], criteria['threshold'], clean_column, clean_column)
            data = data.query(expression)
        for column in self.int_filter_columns:
            clean_column = column.replace(' ', '_')
            data[clean_column] = data[clean_column].astype(str).apply(trim_fraction)
        data.to_csv(self.output_file, sep='\t', header=header, index=False, na_rep='NA')

def trim_fraction(text):
    if '.0' in text:
        return text[:text.rfind('.0')]
    elif 'nan' in text:
        return 'NA'
    return text
