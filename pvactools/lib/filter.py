import pandas as pd

from pvactools.lib.prediction_class import PredictionClass

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
                dtypes = self.set_column_types()
                data = (pd.read_csv(self.input_file, delimiter='\t', float_precision='round_trip', low_memory=False, dtype=dtypes)
                        [lambda x: (x[self.split_column] == key)])
                (filtered_data, header) = self.filter_data(data)
                if append:
                    filtered_data.to_csv(self.output_file, mode='a', sep='\t', header=False, index=False, na_rep='NA')
                else:
                    filtered_data.to_csv(self.output_file, sep='\t', header=header, index=False, na_rep='NA')
                    append = True
        else:
            data = pd.read_csv(self.input_file, delimiter='\t', float_precision='round_trip', low_memory=False)
            (filtered_data, header) = self.filter_data(data)
            filtered_data.to_csv(self.output_file, sep='\t', header=header, index=False, na_rep='NA')

    def get_indexes(self):
        key_df = pd.read_csv(self.input_file, delimiter="\t", usecols=[self.split_column], dtype={self.split_column: str})
        keys = key_df[self.split_column].values.tolist()
        return sorted(list(set(keys)))

    def set_column_types(self):
        dtypes = { self.split_column: str }
        if self.split_column == 'Index':
            dtypes.update({
                'Chromosome': str,
                "Start": "int32",
                "Stop": "int32",
                'Reference': str,
                'Variant': str,
                "Variant Type": "category",
                "Mutation Position": "category",
                "Median MT Score": "float",
                "Median MT Percentile": "float",
                "Protein Position": "str",
            })
            for algorithm in self.determine_used_prediction_algorithms():
                if algorithm == 'SMM' or algorithm == 'SMMPMBEC':
                    continue
                dtypes["{} MT Score".format(algorithm)] = "float"
                dtypes["{} MT Percentile".format(algorithm)] = "float"
        return dtypes

    def determine_used_prediction_algorithms(self):
        headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
        potential_algorithms = PredictionClass.prediction_methods()
        prediction_algorithms = []
        for algorithm in potential_algorithms:
            if "{} MT Score".format(algorithm) in headers or "{} Score".format(algorithm) in headers:
                prediction_algorithms.append(algorithm)
        return prediction_algorithms

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
