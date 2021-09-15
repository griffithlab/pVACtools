import csv
import argparse
import re

import pvactools.lib.sort

class TopScoreFilter:
    def __init__(self, input_file, output_file, top_score_metric, file_type='pVACseq'):
        self.input_file = input_file
        self.output_file = output_file
        self.top_score_metric = top_score_metric
        self.file_type = file_type

    def parse_transcript_id(self, transcript):
        if transcript.startswith('ENS'):
            return re.compile(r'ENS(?:[A-Z]{3})?T(\d+)').match(transcript).group(1)
        elif transcript.startswith('NM_'):
            return re.compile(r'NM_(\d+)').match(transcript).group(0)
        else:
            raise Exception("Error parsing transcript {}. Transcript name doesn't start with ENS or NM".format(transcript))
            return None

    def find_line_with_lowest_transcript_id(self, lines):
        line_with_lowest_transcript_id = lines[0]
        lowest_transcript_id = self.parse_transcript_id(line_with_lowest_transcript_id['Transcript'])
        if lowest_transcript_id is not None:
            for line in lines:
                transcript_id = self.parse_transcript_id(line['Transcript'])
                if transcript_id < lowest_transcript_id:
                    lowest_transcript_id = transcript_id
                    line_with_lowest_transcript_id = line
        return line_with_lowest_transcript_id

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter = "\t")
            writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames = reader.fieldnames)
            writer.writeheader()
            top_per_variant_transcript = {}
            for line in reader:
                if self.file_type == 'pVACseq':
                    chromosome = line['Chromosome']
                    start = line['Start']
                    stop = line['Stop']
                    ref = line['Reference']
                    var = line['Variant']
                    transcript = line['Transcript']
                    index = '%s.%s.%s.%s.%s.%s' % (chromosome, start, stop, ref, var, transcript)
                    if index not in top_per_variant_transcript:
                        top_per_variant_transcript[index] = line
                    top_median_score = float(top_per_variant_transcript[index]['Median MT Score'])
                    top_best_score = float(top_per_variant_transcript[index]['Best MT Score'])
                    median_score = float(line['Median MT Score'])
                    best_score = float(line['Best MT Score'])
                else:
                    index = line['Mutation']
                    if index not in top_per_variant_transcript:
                        top_per_variant_transcript[index] = line
                    top_median_score = float(top_per_variant_transcript[index]['Median Score'])
                    top_best_score = float(top_per_variant_transcript[index]['Best Score'])
                    median_score = float(line['Median Score'])
                    best_score = float(line['Best Score'])
                if ((self.top_score_metric == 'median' and median_score < top_median_score) or
                    (self.top_score_metric == 'lowest' and best_score < top_best_score)):
                    top_per_variant_transcript[index] = line

            top_per_variant = {}
            for (index, line) in top_per_variant_transcript.items():
                if self.file_type == 'pVACseq':
                    chromosome = line['Chromosome']
                    start = line['Start']
                    stop = line['Stop']
                    ref = line['Reference']
                    var = line['Variant']
                    index = '%s.%s.%s.%s.%s' % (chromosome, start, stop, ref, var)
                    epitope = line['MT Epitope Seq']
                else:
                    index = line['Mutation']
                    epitope = line['Epitope Seq']
                if index not in top_per_variant:
                    top_per_variant[index] = {epitope:  [line]}
                else:
                    if epitope in top_per_variant[index]:
                        top_per_variant[index][epitope].append(line)
                    else:
                        top_per_variant[index][epitope] = [line]

            filtered_lines = []
            for (index, per_epitope_lines) in top_per_variant.items():
                for (epitope, lines) in per_epitope_lines.items():
                    if len(lines) == 1:
                        filtered_lines.append(lines[0])
                    else:
                        lines_with_transcript_expression = list(filter(lambda line: line['Transcript Expression'] != 'NA', lines))
                        if len(lines_with_transcript_expression) > 0:
                            line_with_max_expression = lines_with_transcript_expression[0]
                            for line_with_transcript_expression in lines_with_transcript_expression:
                                if float(line_with_transcript_expression['Transcript Expression']) > float(line_with_max_expression['Transcript Expression']):
                                    line_with_max_expression = line_with_transcript_expression
                            filtered_lines.append(line_with_max_expression)
                        else:
                            line_with_lowest_transcript_id = self.find_line_with_lowest_transcript_id(lines)
                            filtered_lines.append(line_with_lowest_transcript_id)

            if self.file_type == 'pVACseq':
                sorted_rows = pvactools.lib.sort.default_sort(filtered_lines, self.top_score_metric)
            else:
                sorted_rows = pvactools.lib.sort.pvacbind_sort(filtered_lines, self.top_score_metric)
            writer.writerows(sorted_rows)

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s top_score_filter' % tool,
            description="Pick the best neoepitope for each variant",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="The final report .tsv file to filter."
        )
        parser.add_argument(
            'output_file',
            help="Output .tsv file containing only the list of the top "
                 + "epitope per variant."
        )
        parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            default='median',
            help="The ic50 scoring metric to use for filtering. "
                 + "lowest: Use the best MT Score (i.e. the lowest MT ic50 binding score of all chosen prediction methods). "
                 + "median: Use the median MT Score (i.e. the median MT ic50 binding score of all chosen prediction methods)."
        )
        return parser
