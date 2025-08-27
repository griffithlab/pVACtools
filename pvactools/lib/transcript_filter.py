import argparse
import sys
import re
import csv
from pvactools.lib.run_utils import *

class TranscriptFilter:
    def __init__(self, input_file, output_file, transcript_prioritization_strategy= ['canonical', 'mane_select', 'tsl'], maximum_transcript_support_level=1):
        self.input_file = input_file
        self.output_file = output_file
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter="\t")
            writer = csv.DictWriter(output_fh, delimiter="\t", fieldnames=reader.fieldnames)
            writer.writeheader()
            output_lines = []
            for line in reader:
                transcript_pass = is_preferred_transcript(line, self.transcript_prioritization_strategy, self.maximum_transcript_support_level)
                if transcript_pass:
                    output_lines.append(line)
            writer.writerows(output_lines)

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s transcript_filter' % tool,
            description="Filter variant transcripts processed by IEDB.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="The all_epitopes.tsv or filtered.tsv report file to filter."
        )
        parser.add_argument(
            'output_file',
            help="Output .tsv file containing list of filtered "
                 + "epitopes based on the variant transcript."
        )
        parser.add_argument(
            "--transcript-prioritization-strategy", type=transcript_prioritization_strategy(),
            help="Specify the criteria to consider when filtering transcripts of the neoantigen candidates. "
                 + "'canonical' will select candidates resulting from variants on a Ensembl canonical transcript. "
                 + "'mane_select' will select candidates resulting from variants on a MANE select transcript. "
                 + "'tsl' will select candidates where the transcript support level (TSL) matches the --maximum-transcript-support-level cutoff. "
                 + "When selecting more than one criteria, a transcript meeting EITHER of the selected criteria will be selected.",
            default=['canonical', 'mane_select', 'tsl']
        )
        parser.add_argument(
            "--maximum-transcript-support-level", type=int,
            help="The threshold to use for filtering epitopes on the Ensembl transcript support level (TSL). "
                 + "Keep all epitopes with a transcript support level <= to this cutoff.",
            default=1,
            choices=[1, 2, 3, 4, 5]
        )
        return parser
