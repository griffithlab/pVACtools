import argparse
import sys
import re
import csv
import json

from pvactools.lib.run_argument_utils import tiers

class AggregateReportFilter:
    def __init__(self, input_file, output_file, input_metrics_file=None, output_metrics_file=None, include_tiers=["Pass"]):
        self.input_file = input_file
        self.output_file = output_file
        self.include_tiers = include_tiers
        self.input_metrics_file = input_metrics_file
        self.output_metrics_file = output_metrics_file

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter="\t")
            writer = csv.DictWriter(output_fh, delimiter="\t", fieldnames=reader.fieldnames)
            writer.writeheader()
            output_lines = []
            remove_lines = []
            for line in reader:
                if line['Tier'] in self.include_tiers:
                    output_lines.append(line)
                else:
                    remove_lines.append(line)
            if self.input_metrics_file and self.output_metrics_file:
                with open(self.input_metrics_file, 'r') as fh:
                    aggregate_metrics = json.loads(fh.read())
                remove_keys = [line['ID'] for line in remove_lines]
                for key in remove_keys:
                    del aggregate_metrics[key]
                with open(self.output_metrics_file, 'w') as fh:
                    json.dump(aggregate_metrics, fh, indent=2, separators=(',', ': '))
            writer.writerows(output_lines)

    @classmethod
    def parser(cls, tool):
        if tool == 'pvacseq':
            description="Filter an aggregate report and its metrics.json file based on the variant Tier."
        else:
            description="Filter an aggregate report based on the variant Tier."
        parser = argparse.ArgumentParser(
            '%s aggregate_report_filter' % tool,
            description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="The aggregated.tsv report file to filter."
        )
        parser.add_argument(
            'output_file',
            help="Output aggregated.tsv file containing list of filtered "
                 + "aggregate report entries based on the selected variant Tier."
        )
        if tool == 'pvacseq':
            parser.add_argument(
                'input_metrics_file',
                help="The metrics.json file accompanying the input aggregated report file."
            )
            parser.add_argument(
                'output_metrics_file',
                help="Filtered metrics.json file only retaining those entries from variants in the selected variant Tier."
            )
        if tool == 'pvacseq':
            help_text = "Specify a comma-separated list of tiers for which to retain aggregate report and metrics.json variant data."
        else:
            help_text = "Specify a comma-separated list of tiers for which to retain aggregate report variant data."
        parser.add_argument(
            "--include-tiers", type=tiers(tool),
            help=help_text,
            default=['Pass']
        )
        return parser
