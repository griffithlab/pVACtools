import argparse
import sys
import re
import csv
import textwrap
import pandas as pd

from pvactools.lib.run_utils import *


class IdentifyProblematicAminoAcids:
    def __init__(self, input_file, output_file, problematic_amino_acids, file_type='pVACseq', filter_type='soft'):
        self.input_file = input_file
        self.output_file = output_file
        self.problematic_amino_acids = [IdentifyProblematicAminoAcids.process_problematic_aa_entry(e) for e in problematic_amino_acids]
        self.file_type = file_type
        self.filter_type = filter_type
        self.is_aggregated = self._is_aggregated()
        self.prob_column_exists = self._prob_column_exists()

    def _is_aggregated(self):
        headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
        return 'Best Peptide' in headers

    def _prob_column_exists(self):
        headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
        if self.is_aggregated:
            return 'Prob Pos' in headers
        else:
            return 'Problematic Positions' in headers

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter = "\t")
            if self.filter_type == 'soft':
                if self.prob_column_exists:
                    writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=reader.fieldnames, extrasaction='ignore', restval='NA')
                else:
                    if self.is_aggregated:
                        writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=reader.fieldnames + ["Prob Pos"], extrasaction='ignore', restval='NA')
                    else:
                        writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=reader.fieldnames + ["Problematic Positions"], extrasaction='ignore', restval='NA')
            else:
                writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=reader.fieldnames, extrasaction='ignore', restval='NA')
            writer.writeheader()
            for line in reader:
                if self.is_aggregated:
                    sequence = line['Best Peptide']
                    output_header = 'Prob Pos'
                elif self.file_type == 'pVACbind' or self.file_type == 'pVACfuse' or self.file_type == 'pVACsplice':
                    sequence = line['Epitope Seq']
                    output_header = 'Problematic Positions'
                else:
                    sequence = line['MT Epitope Seq']
                    output_header = 'Problematic Positions'
                problematic_positions = []
                for (aa, position) in self.problematic_amino_acids:
                    if position is None and aa in sequence:
                        problematic_positions.extend([str(m.start()+1) for m in re.finditer(aa, sequence)])
                    elif position is not None and aa in sequence and sequence.index(aa) + 1 == position:
                        problematic_positions.append(str(position))
                if self.filter_type == 'hard' and len(problematic_positions) == 0:
                    writer.writerow(line)
                elif self.filter_type == 'soft':
                    if len(problematic_positions) == 0:
                        line[output_header] = 'None'
                    else:
                        line[output_header] = ','.join(problematic_positions)
                    writer.writerow(line)

    @classmethod
    def process_problematic_aa_entry(cls, problematic_amino_acid):
        supported_aas = supported_amino_acids()
        if ":" in problematic_amino_acid:
            (amino_acids, position) = problematic_amino_acid.split(":")
            for aa in amino_acids:
                if not aa in supported_aas:
                    raise Exception("Error with problematic amino acid input entry \"{}\". Amino acid {} not supported.".format(problematic_amino_acid, aa))
            try:
                position = int(position)
            except:
                raise Exception("Error with problematic amino acid input entry \"{}\". Position {} isn't an integer.".format(problematic_amino_acid, position))
            if position == 0:
                raise Exception("Error with problematic amino acid input entry \"{}\". Position can't be 0.".format(problematic_amino_acid, position))
            return (amino_acids, position)
        else:
            for aa in problematic_amino_acid:
                if not aa in supported_aas:
                    raise Exception("Error with problematic amino acid input entry \"{}\". Amino acid {} not supported.".format(problematic_amino_acid, aa))
            return (problematic_amino_acid, None)

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s identify_problematic_amino_acids' % tool,
            description="Mark problematic amino acid positions in each epitope or filter entries that have problematic amino acids.",
            formatter_class=argparse.RawTextHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="Input filtered, all_epitopes, or aggregated file with predicted epitopes."
        )
        parser.add_argument(
            'output_file',
            help="Output .tsv file with identification of problematic amino acids or hard-filtered to remove epitopes with problematic amino acids."
        )
        parser.add_argument(
            'problematic_amino_acids', type=lambda s:[a for a in s.split(',')],
            help=textwrap.dedent('''\
            A list of amino acids to consider as problematic. Each entry can be specified in the following format:
            `amino_acid(s)`: One or more one-letter amino acid codes. Any occurrence of this amino acid string,
                             regardless of the position in the epitope, is problematic. When specifying more than
                             one amino acid, they will need to occur together in the specified order.
            `amino_acid:position`: A one letter amino acid code, followed by a colon separator, followed by a positive
                                   integer position (one-based). The occurrence of this amino acid at the position
                                   specified is problematic., E.g. G:2 would check for a Glycine at the second position
                                   of the epitope. The N-terminus is defined as position 1.
            `amino_acid:-position`: A one letter amino acid code, followed by a colon separator, followed by a negative
                                    integer position. The occurrence of this amino acid at the specified position from
                                    the end of the epitope is problematic. E.g., G:-3 would check for a Glycine at the
                                    third position from the end of the epitope. The C-terminus is defined as position -1.''')
        )
        parser.add_argument(
            '--filter-type', '-f', choices=['soft', 'hard'], default="soft",
            help="Set the type of filtering done. Choosing `soft` will add a new column \"Problematic Positions\" (for filtered or all_epitopes input files) "
               + "or \"Prob Pos\" (for aggregated input files) that lists positions in the epitope with problematic amino acids. "
               + "Choosing `hard` will remove epitope entries with problematic amino acids."
        )
        return parser
