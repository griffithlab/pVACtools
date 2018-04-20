import csv

class CondenseFinalReport:
    def __init__(self, input_file, output_file, top_score_metric):
        self.input_file = input_file
        self.output_file = output_file
        self.top_score_metric = top_score_metric

    def condensed_header(self):
        return [
            'Gene Name',
            'Mutation',
            'Protein Position',
            'HGVSc',
            'HGVSp',
            'HLA Allele',
            'MT Epitope Seq',
            'MT IC50',
            'WT IC50',
            'Fold Change',
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
                if self.top_score_metric == 'median':
                    line['MT IC50'] = line['Median MT Score']
                    line['WT IC50'] = line['Median WT Score']
                    line['Fold Change'] = line['Median Fold Change']
                elif self.top_score_metric == 'lowest':
                    line['MT IC50'] = line['Best MT Score']
                    line['WT IC50'] = line['Corresponding WT Score']
                    line['Fold Change'] = line['Corresponding Fold Change']
                writer.writerow(line)
