import argparse
import os
import csv


class MarkGenesOfInterest:
    def __init__(self, input_file, output_file, genes_of_interest_file=None):
        self.input_file = input_file
        self.output_file = output_file
        if genes_of_interest_file is None:
            base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
            genes_of_interest_file = os.path.join(base_dir, 'supporting_files', 'cancer_census_hotspot_gene_list.tsv')

        with open(genes_of_interest_file) as fh:
            self.genes_of_interest = [line.rstrip() for line in fh]

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter = "\t")
            writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=reader.fieldnames + ["Gene of Interest"], extrasaction='ignore', restval='NA')
            writer.writeheader()
            for line in reader:
                if line['Gene Name'] in self.genes_of_interest:
                    line['Gene of Interest'] = "True"
                else:
                    line['Gene of Interest'] = "False"
                writer.writerow(line)

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s mark_genes_of_interest' % tool,
            description="Mark predictions resulting from variants on a genes of interest list.",
            formatter_class=argparse.RawTextHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="Input filtered or all_epitopes file with predicted epitopes."
        )
        parser.add_argument(
            'output_file',
            help="Output .tsv file with identification of predictions resulting from variants on the genes of interest list."
        )
        parser.add_argument(
            '--genes-of-interest-file',
            help="A genes of interest file. Predictions resulting from variants on genes in this list will be marked in the output file. "
                 + "The file should be formatted to have each gene on a separate line without a header line. "
                 + "If no file is specified, the Cancer Gene Census list of high-confidence genes is used as the default."
        )
        return parser
