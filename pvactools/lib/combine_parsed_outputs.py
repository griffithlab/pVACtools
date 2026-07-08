import sys
import csv
import pvactools.lib.sort

class CombineParsedOutputs:
    def __init__(self, **kwargs):
        self.input_files = kwargs['input_files']
        self.output_file = kwargs['output_file']

    def execute(self):
        fieldnames = []
        for input_file in self.input_files:
            with open(input_file, 'r') as input_file_handle:
                reader = csv.DictReader(input_file_handle, delimiter='\t')
                if len(fieldnames) == 0:
                    fieldnames = reader.fieldnames
                else:
                    for fieldname in reader.fieldnames:
                        if fieldname not in fieldnames:
                            fieldnames.append(fieldname)

        with open(self.output_file, 'w') as fh:
            tsv_writer = csv.DictWriter(fh, list(fieldnames), delimiter = '\t', lineterminator = '\n', restval='NA')
            tsv_writer.writeheader()
            for input_file in self.input_files:
                with open(input_file, 'r') as input_file_handle:
                    reader = csv.DictReader(input_file_handle, delimiter='\t')
                    for row in reader:
                        tsv_writer.writerow(row)
