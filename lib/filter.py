import pandas

class Filter:
    def __init__(self, input_file, output_file, filter_criteria, exclude_nas):
        self.input_file = input_file
        self.output_file = output_file
        self.filter_criteria = filter_criteria
        self.exclude_nas = exclude_nas

    def execute(self):
        data = pandas.read_csv(self.input_file, delimiter='\t', float_precision='high')
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
        data.to_csv(self.output_file, sep='\t', header=header, index=False, na_rep='NA')
