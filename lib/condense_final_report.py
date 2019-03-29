import csv

class CondenseFinalReport:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def condensed_header(self):
        return [
            'Gene Name',
            'Mutation',
            'Protein Position',
            'HGVSc',
            'HGVSp',
            'HLA Allele',
            'Mutation Position',
            'MT Epitope Seq',
            'Median MT Score',
            'Median WT Score',
            'Median Fold Change',
            'Best MT Score',
            'Corresponding WT Score',
            'Corresponding Fold Change',
            'Tumor DNA Depth',
            'Tumor DNA VAF',
            'Tumor RNA Depth',
            'Tumor RNA VAF',
            'Gene Expression',
        ]

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter = "\t")
            writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=self.condensed_header(), extrasaction='ignore')
            writer.writeheader()

            for line in reader:
                writer.writerow(line)
