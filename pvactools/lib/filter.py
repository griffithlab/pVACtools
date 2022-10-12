import pandas as pd

pd.options.mode.chained_assignment = None

class Filter:
    def __init__(self, input_file, output_file, filter_criteria, int_filter_columns=[], split=False, split_column=None):
        self.input_file = input_file
        self.output_file = output_file
        self.filter_criteria = filter_criteria
        self.int_filter_columns = int_filter_columns
        self.split = split
        self.split_column = split_column

    def execute(self):
        if self.split:
            keys = self.get_indexes()
            append = False
            for key in keys:
                data = (pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False, dtype={self.split_column: str})
                        [lambda x: (x[self.split_column] == key)])
                (filtered_data, header) = self.filter_data(data)
                if append:
                    filtered_data.to_csv(self.output_file, mode='a', sep='\t', header=False, index=False, na_rep='NA')
                else:
                    filtered_data.to_csv(self.output_file, sep='\t', header=header, index=False, na_rep='NA')
                    append = True
        else:
            data = pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False)
            (filtered_data, header) = self.filter_data(data)
            filtered_data.to_csv(self.output_file, sep='\t', header=header, index=False, na_rep='NA')

    def get_indexes(self):
        key_df = pd.read_csv(self.input_file, delimiter="\t", usecols=[self.split_column], dtype={self.split_column: str})
        keys = key_df[self.split_column].values.tolist()
        return sorted(list(set(keys)))

    def filter_data(self, data):
        header = data.columns
        clean_header = header.map(lambda x: x.replace(' ', '_') if isinstance(x, str) else x)
        data.columns = clean_header
        for criteria in self.filter_criteria:
            clean_column = criteria['column'].replace(' ', '_')
            #clean_column != clean_column is a hacky way to keep all NA values
            if criteria['exclude_nas']:
                expression = "(%s %s %s)" % (clean_column, criteria['operator'], criteria['threshold'])
            else:
                expression = "(%s %s %s) | (%s != %s)" % (clean_column, criteria['operator'], criteria['threshold'], clean_column, clean_column)
            data = data.query(expression)
        for column in self.int_filter_columns:
            clean_column = column.replace(' ', '_')
            data[clean_column] = data[clean_column].astype(str).apply(trim_fraction)
        return (data, header)

def trim_fraction(text):
    if '.0' in text:
        return text[:text.rfind('.0')]
    elif 'nan' in text:
        return 'NA'
    return text
