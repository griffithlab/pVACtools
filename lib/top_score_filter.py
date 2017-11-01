import csv
import argparse

class TopScoreFilter:
    def __init__(self, input_file, output_file, top_score_metric):
        self.input_file = input_file
        self.output_file = output_file
        self.top_score_metric = top_score_metric

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter = "\t")
            writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames = reader.fieldnames)
            writer.writeheader()
            filtered_results = {}
            for line in reader:
                gene_name = line['Gene Name']
                transcript_name = line['Transcript']
                consequence = line['Variant Type']
                if consequence == 'FS':
                    amino_acid_change_position = "%s%s/%s" % (line['Protein Position'], line['Reference'], line ['Variant'])
                elif 'fusion' in consequence:
                    amino_acid_change_position = line['Fusion Position']
                else:
                    amino_acid_change_position = line['Protein Position'] + line['Mutation']
                index = '%s.%s.%s.%s' % (gene_name, transcript_name, consequence, amino_acid_change_position)
                if index not in filtered_results:
                    filtered_results[index] = line
                else:
                    if ((self.top_score_metric == 'median' and float(line['Median MT Score']) < float(filtered_results[index]['Median MT Score'])) or
                        (self.top_score_metric == 'lowest' and float(line['Best MT Score']) < float(filtered_results[index]['Best MT Score']))):
                        filtered_results[index] = line

            writer.writerows(filtered_results.values())

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser('%s top_score_filter' % tool)
        parser.add_argument(
            'input_file',
            help="The final report .tsv file to filter"
        )
        parser.add_argument(
            'output_file',
            help="Output .tsv file containing only the list of the top "
                 + "epitope per variant"
        )
        parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            default='median',
            help="The ic50 scoring metric to use for filtering. "
                 + "lowest: Best MT Score - lowest MT ic50 binding score of all chosen prediction methods. "
                 + "median: Median MT Score - median MT ic50 binding score of all chosen prediction methods. "
                 + "Default: median"
        )
        return parser
