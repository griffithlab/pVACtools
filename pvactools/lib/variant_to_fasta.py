import csv
import json
import re
import sys
from collections import OrderedDict, defaultdict
import yaml
from abc import ABCMeta
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import itertools
import logging

from pvactools.lib.proximal_variant import ProximalVariant
from pvactools.lib.run_utils import *

csv.field_size_limit(sys.maxsize)

class VariantToFasta(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        self.input_file                 = kwargs['input_file']
        self.output_file                = kwargs['output_file']
        self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', None)
        self.proximal_variants_file     = kwargs.pop('proximal_variants_file', None)
        self.proximal_variants          = self.parse_proximal_variants_file()

    def parse_proximal_variants_file(self):
        if self.proximal_variants_file is not None:
            proximal_variants = defaultdict(lambda: defaultdict(list))
            with open(self.proximal_variants_file, 'r') as fh:
                tsvin = csv.DictReader(fh, delimiter='\t')
                for line in tsvin:
                    proximal_variants[line['main_somatic_variant']][line['protein_position']].append(line)
            return proximal_variants
        else:
            return {}

    def invalid_characters(self):
        return ['*', 'X', '?']

    def position_out_of_bounds(self, position, sequence):
        return position > len(sequence)-1

    def add_proximal_variants(self, somatic_variant_index, wildtype_sequence, position, germline_variants_only):
        wildtype_sequence_with_proximal_variants = wildtype_sequence
        if somatic_variant_index in self.proximal_variants.keys():
            for (protein_position, lines) in self.proximal_variants[somatic_variant_index].items():
                if protein_position == position:
                    continue

                if germline_variants_only:
                    filtered_lines = [line for line in lines if line['type'] == 'germline']
                else:
                    filtered_lines = lines

                if len(filtered_lines) == 0:
                    continue
                elif len(filtered_lines) == 1:
                    line = filtered_lines[0]
                    proximal_variant_wildtype_amino_acid, proximal_variant_mutant_amino_acid = line['amino_acid_change'].split('/')
                else:
                    line = filtered_lines[0]
                    proximal_variant_wildtype_amino_acid = line['amino_acid_change'].split('/')[0]
                    codon_changes = [ item['codon_change'] for item in filtered_lines ]
                    proximal_variant_mutant_amino_acid = ProximalVariant.combine_conflicting_variants(codon_changes)

                if '-' in protein_position:
                    proximal_variant_start_position = int(protein_position.split('-')[0]) - 1
                else:
                    proximal_variant_start_position = int(protein_position) - 1
                proximal_variant_end_position = proximal_variant_start_position + len(proximal_variant_mutant_amino_acid)

                if proximal_variant_end_position <= 0 or proximal_variant_start_position >= len(wildtype_sequence):
                    continue
                if len(proximal_variant_wildtype_amino_acid) != len(proximal_variant_mutant_amino_acid):
                    print("Nearby variant is not a missense mutation. Skipping.")
                    continue

                if proximal_variant_end_position > len(wildtype_sequence):
                    #the DNP extends past the end of the wildtype_sequence
                    missing_amino_acids_count = proximal_variant_end_position - len(wildtype_sequence)
                    #remove proximal variant amino acids after wildtype subsquence end
                    proximal_variant_wildtype_amino_acid = proximal_variant_wildtype_amino_acid[:len(proximal_variant_wildtype_amino_acid) - missing_amino_acids_count]
                    proximal_variant_mutant_amino_acid = proximal_variant_mutant_amino_acid[:len(proximal_variant_mutant_amino_acid) - missing_amino_acids_count]

                if proximal_variant_start_position < 0:
                    #DNP starts before the beginning of the wildtype sequence
                    missing_amino_acids_count = abs(proximal_variant_start_position)
                    proximal_variant_start_position = 0
                    #remove proximal variant amino acids before wildtype sequence start
                    proximal_variant_wildtype_amino_acid = proximal_variant_wildtype_amino_acid[missing_amino_acids_count:]
                    proximal_variant_mutant_amino_acid = proximal_variant_mutant_amino_acid[missing_amino_acids_count:]

                if wildtype_sequence[proximal_variant_start_position:proximal_variant_end_position] != proximal_variant_wildtype_amino_acid:
                    sys.exit(
                        "Error when processing proximal variant.\n" +
                        "The wildtype amino acid for variant %s with substring %s is different than expected.\n" % (somatic_variant_index, wildtype_sequence) +
                        "Actual wildtype amino acid: %s\n" % wildtype_sequence[proximal_variant_start_position:proximal_variant_end_position] +
                        "Wildtype amino acid of the proximal_variant: %s" % proximal_variant_wildtype_amino_acid
                    )

                wildtype_sequence_with_proximal_variants = wildtype_sequence_with_proximal_variants[:proximal_variant_start_position] + proximal_variant_mutant_amino_acid + wildtype_sequence_with_proximal_variants[proximal_variant_end_position:]
        return wildtype_sequence_with_proximal_variants

    def execute(self):
        records = []
        with open(self.input_file, 'r') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for line in reader:
                variant_type = line['variant_type']
                wildtype_sequence = line['wildtype_amino_acid_sequence']
                if variant_type == 'FS':
                    position = int(line['protein_position'].split('-', 1)[0]) - 1
                elif variant_type == 'missense' or variant_type == 'inframe_ins':
                    if '/' not in line['amino_acid_change']:
                        continue
                    wildtype_amino_acid, mutant_amino_acid = line['amino_acid_change'].split('/')
                    if '*' in wildtype_amino_acid:
                        wildtype_amino_acid = wildtype_amino_acid.split('*')[0]
                    elif 'X' in wildtype_amino_acid:
                        wildtype_amino_acid = wildtype_amino_acid.split('X')[0]
                    if '*' in mutant_amino_acid:
                        mutant_amino_acid = mutant_amino_acid.split('*')[0]
                        stop_codon_added = True
                    elif 'X' in mutant_amino_acid:
                        mutant_amino_acid = mutant_amino_acid.split('X')[0]
                        stop_codon_added = True
                    else:
                        stop_codon_added = False
                    if wildtype_amino_acid == '-':
                        position = int(line['protein_position'].split('-', 1)[0])
                        wildtype_amino_acid_length = 0
                    else:
                        if '-' in line['protein_position']:
                            position = int(line['protein_position'].split('-', 1)[0]) - 1
                            wildtype_amino_acid_length = len(wildtype_amino_acid)
                        else:
                            position = int(line['protein_position']) - 1
                            wildtype_amino_acid_length = len(wildtype_amino_acid)
                elif variant_type == 'inframe_del':
                    variant_type = 'inframe_del'
                    wildtype_amino_acid, mutant_amino_acid = line['amino_acid_change'].split('/')
                    if '*' in wildtype_amino_acid:
                        wildtype_amino_acid = wildtype_amino_acid.split('*')[0]
                    elif 'X' in wildtype_amino_acid:
                        wildtype_amino_acid = wildtype_amino_acid.split('X')[0]
                    if '*' in mutant_amino_acid:
                        mutant_amino_acid = mutant_amino_acid.split('*')[0]
                        stop_codon_added = True
                    elif 'X' in mutant_amino_acid:
                        mutant_amino_acid = mutant_amino_acid.split('X')[0]
                        stop_codon_added = True
                    else:
                        stop_codon_added = False
                    position = int(line['protein_position'].split('-', 1)[0]) - 1
                    wildtype_amino_acid_length = len(wildtype_amino_acid)
                    if mutant_amino_acid == '-':
                        mutant_amino_acid = ''
                else:
                    continue

                if self.position_out_of_bounds(position, wildtype_sequence):
                    continue

                if variant_type == 'missense' and line['index'] in self.proximal_variants and line['protein_position'] in self.proximal_variants[line['index']]:
                    codon_changes = [ item['codon_change'] for item in self.proximal_variants[line['index']][line['protein_position']] ]
                    codon_changes.append(line['codon_change'])
                    mutant_amino_acid_with_proximal_variants = ProximalVariant.combine_conflicting_variants(codon_changes)
                elif variant_type != 'FS':
                    mutant_amino_acid_with_proximal_variants = mutant_amino_acid

                if variant_type == 'FS':
                    mutant_sequence = line['frameshift_amino_acid_sequence']
                    if self.downstream_sequence_length is not None:
                        mutant_sequence = mutant_sequence[:(position + self.downstream_sequence_length)]
                    if wildtype_sequence.startswith(mutant_sequence):
                        continue
                    left_flanking_sequence = mutant_sequence[:position]
                    wildtype_sequence = self.add_proximal_variants(line['index'], wildtype_sequence, position, True)
                    left_flanking_sequence_with_proximal_variants = self.add_proximal_variants(line['index'], left_flanking_sequence, position, False)
                    #The caveat here is that if a nearby variant is in the downstream sequence, the protein sequence would be further altered, which we aren't taking into account.
                    #we would need to recalculate the downstream protein sequence taking all downstream variants into account.
                    mutant_sequence = re.sub('^%s' % left_flanking_sequence, left_flanking_sequence_with_proximal_variants, mutant_sequence)
                else:
                    end_position = position + wildtype_amino_acid_length
                    if wildtype_amino_acid != '-' and wildtype_amino_acid != wildtype_sequence[position:end_position]:
                        if line['amino_acid_change'].split('/')[0].count('*') > 1:
                            print("Warning: Amino acid change is not sane - contains multiple stops. Skipping entry {}".format(line['index']))
                            continue
                        else:
                            sys.exit("ERROR: There was a mismatch between the actual wildtype amino acid sequence ({}) and the expected amino acid sequence ({}). Did you use the same reference build version for VEP that you used for creating the VCF?\n{}".format(wildtype_sequence[position:end_position], wildtype_amino_acid, line))
                    wildtype_sequence_with_proximal_variants = self.add_proximal_variants(line['index'], wildtype_sequence, position, False)
                    wildtype_sequence = self.add_proximal_variants(line['index'], wildtype_sequence, position, True)
                    if stop_codon_added:
                        mutant_sequence = wildtype_sequence_with_proximal_variants[:position] + mutant_amino_acid_with_proximal_variants
                    else:
                        mutant_sequence = wildtype_sequence_with_proximal_variants[:position] + mutant_amino_acid_with_proximal_variants + wildtype_sequence_with_proximal_variants[end_position:]

                if mutant_sequence in wildtype_sequence:
                    continue

                index = line['index']
                records.append(SeqRecord(Seq(mutant_sequence), id=f"MT.{index}", description=""))
                records.append(SeqRecord(Seq(wildtype_sequence), id=f"WT.{index}", description=""))

            SeqIO.write(records, self.output_file, "fasta")
