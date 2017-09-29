import pandas

class Filter:
    def __init__(self, input_file, output_file, filter_criteria):
        self.input_file = input_file
        self.output_file = output_file
        self.filter_criteria = filter_criteria

    def execute(self):
        data = pandas.read_csv(self.input_file, delimiter='\t', float_precision='high')
        header = data.columns
        clean_header = header.map(lambda x: x.replace(' ', '_') if isinstance(x, str) else x)
        data.columns = clean_header
        for criteria in self.filter_criteria:
            data = data.query(criteria)
        data.to_csv(self.output_file, sep='\t', header=header, index=False, na_rep='NA')
