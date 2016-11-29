import csv
import re
import sys
from collections import OrderedDict
import yaml
from abc import ABCMeta

csv.field_size_limit(sys.maxsize)

class FastaGenerator(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        self.input_file                 = kwargs['input_file']
        self.peptide_sequence_length    = kwargs['peptide_sequence_length']
        self.epitope_length             = kwargs['epitope_length']
        self.output_file                = kwargs['output_file']
        self.output_key_file            = kwargs['output_key_file']
        self.downstream_sequence_length = kwargs['downstream_sequence_length']

    def position_out_of_bounds(self, position, sequence):
        return position > len(sequence)-1

    #This subroutine is a bit funky but it was designed that way to mirror
    #distance_from_end to increase code readability from the caller's perspective
    def distance_from_start(self, position, string):
        return position

    def distance_from_end(self, position, string):
        return len(string) - 1 - position

    def determine_peptide_sequence_length(self, full_wildtype_sequence_length, peptide_sequence_length, line):
        actual_peptide_sequence_length = peptide_sequence_length

        #If the wildtype sequence is shorter than the desired peptide sequence
        #length we use the wildtype sequence length instead so that the extraction
        #algorithm below works correctly
        if full_wildtype_sequence_length < actual_peptide_sequence_length:
            actual_peptide_sequence_length = full_wildtype_sequence_length
            print('Wildtype sequence length is shorter than desired peptide sequence length at position (%s, %s, %s). Using wildtype sequence length (%s) instead.' % (line['chromosome_name'], line['start'], line['stop'], actual_peptide_sequence_length))

        return actual_peptide_sequence_length

    def determine_flanking_sequence_length(self, full_wildtype_sequence_length, peptide_sequence_length, line):
        actual_peptide_sequence_length = self.determine_peptide_sequence_length(full_wildtype_sequence_length, peptide_sequence_length, line)
        if actual_peptide_sequence_length%2 == 0:
            return int((actual_peptide_sequence_length-2) / 2)
        else:
            return int((actual_peptide_sequence_length-1) / 2)

    def get_wildtype_subsequence(self, position, full_wildtype_sequence, wildtype_amino_acid_length, peptide_sequence_length, line):
        one_flanking_sequence_length = self.determine_flanking_sequence_length(len(full_wildtype_sequence), peptide_sequence_length, line)
        peptide_sequence_length = 2 * one_flanking_sequence_length + wildtype_amino_acid_length

        # We want to extract a subset from full_wildtype_sequence that is
        # peptide_sequence_length long so that the position ends
        # up in the middle of the extracted sequence.
        # If the position is too far toward the beginning or end of
        # full_wildtype_sequence there aren't enough amino acids on one side
        # to achieve this.
        if self.distance_from_start(position, full_wildtype_sequence) < one_flanking_sequence_length:
            wildtype_subsequence = full_wildtype_sequence[:peptide_sequence_length]
            mutation_position = position
        elif self.distance_from_end(position, full_wildtype_sequence) < one_flanking_sequence_length:
            start_position = len(full_wildtype_sequence) - peptide_sequence_length
            wildtype_subsequence = full_wildtype_sequence[start_position:]
            mutation_position = peptide_sequence_length - self.distance_from_end(position, full_wildtype_sequence) - 1
        elif self.distance_from_start(position, full_wildtype_sequence) >= one_flanking_sequence_length and self.distance_from_end(position, full_wildtype_sequence) >= one_flanking_sequence_length:
            start_position = position - one_flanking_sequence_length
            end_position   = start_position + peptide_sequence_length
            wildtype_subsequence = full_wildtype_sequence[start_position:end_position]
            mutation_position = one_flanking_sequence_length
        else:
            sys.exit("ERROR: Something went wrong during the retrieval of the wildtype sequence at position(%s, %s, %s)" % line['chromsome_name'], line['start'], line['stop'])

        return mutation_position, wildtype_subsequence

    def get_frameshift_subsequences(self, position, full_wildtype_sequence, peptide_sequence_length, line):
        one_flanking_sequence_length = self.determine_flanking_sequence_length(len(full_wildtype_sequence), peptide_sequence_length, line)
        if position < one_flanking_sequence_length:
            start_position = 0
        else:
            start_position = position - one_flanking_sequence_length
        wildtype_subsequence_stop_position = position + one_flanking_sequence_length
        mutation_subsequence_stop_position = position
        wildtype_subsequence = full_wildtype_sequence[start_position:wildtype_subsequence_stop_position]
        mutation_start_subsequence = full_wildtype_sequence[start_position:mutation_subsequence_stop_position]
        return wildtype_subsequence, mutation_start_subsequence

    def execute(self):
        peptide_sequence_length = self.peptide_sequence_length
        reader                  = open(self.input_file, 'r')
        tsvin                   = csv.DictReader(reader, delimiter='\t')
        fasta_sequences         = OrderedDict()
        for line in tsvin:
            variant_type = line['variant_type']
            full_wildtype_sequence = line['wildtype_amino_acid_sequence']
            if variant_type == 'FS':
                position = int(line['protein_position'].split('-', 1)[0]) - 1
            elif variant_type == 'missense' or variant_type == 'inframe_ins':
                wildtype_amino_acid, mutant_amino_acid = line['amino_acid_change'].split('/')
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
                position = int(line['protein_position'].split('-', 1)[0]) - 1
                wildtype_amino_acid_length = len(wildtype_amino_acid)
                if mutant_amino_acid == '-':
                    mutant_amino_acid = ''
            else:
                continue

            if self.position_out_of_bounds(position, full_wildtype_sequence):
                continue

            if variant_type == 'FS':
                wildtype_subsequence, mutant_subsequence = self.get_frameshift_subsequences(position, full_wildtype_sequence, peptide_sequence_length, line)
                downstream_sequence = line['downstream_amino_acid_sequence']

                if self.downstream_sequence_length and len(downstream_sequence) > self.downstream_sequence_length:
                    downstream_sequence = downstream_sequence[0:self.downstream_sequence_length]
                mutant_subsequence += downstream_sequence
            else:
                mutation_start_position, wildtype_subsequence = self.get_wildtype_subsequence(position, full_wildtype_sequence, wildtype_amino_acid_length, peptide_sequence_length, line)
                mutation_end_position = mutation_start_position + wildtype_amino_acid_length
                mutant_subsequence = wildtype_subsequence[:mutation_start_position] + mutant_amino_acid + wildtype_subsequence[mutation_end_position:]

            if '*' in wildtype_subsequence or '*' in mutant_subsequence:
                continue

            if 'X' in wildtype_subsequence or 'X' in mutant_subsequence:
                continue

            if len(wildtype_subsequence) < self.epitope_length or len(mutant_subsequence) < self.epitope_length:
                continue

            variant_id = line['index']
            for designation, subsequence in zip(['WT', 'MT'], [wildtype_subsequence, mutant_subsequence]):
                key = '%s.%s' % (designation, variant_id)
                if subsequence in fasta_sequences:
                    fasta_sequences[subsequence].append(key)
                else:
                    fasta_sequences[subsequence] = [key]

        writer = open(self.output_file, 'w')
        key_writer = open(self.output_key_file, 'w')
        count  = 1
        for (subsequence, keys) in fasta_sequences.items():
            writer.writelines('>%s\n' % count)
            writer.writelines('%s\n' % subsequence)
            yaml.dump({count: keys}, key_writer, default_flow_style=False)
            count += 1

        reader.close()
        writer.close()
        key_writer.close()

class FusionFastaGenerator(FastaGenerator):
    def execute(self):
        peptide_sequence_length = self.peptide_sequence_length
        reader                  = open(self.input_file, 'r')
        tsvin                   = csv.DictReader(reader, delimiter='\t')
        fasta_sequences         = OrderedDict()
        for line in tsvin:
            variant_type = line['variant_type']
            position     = int(line['fusion_position'])
            sequence     = line['fusion_amino_acid_sequence']
            one_flanking_sequence_length = self.determine_flanking_sequence_length(len(sequence), peptide_sequence_length, line)
            if position < one_flanking_sequence_length:
                start_position = 0
            else:
                start_position = position - one_flanking_sequence_length

            if variant_type == 'inframe_fusion':
                stop_position = position + one_flanking_sequence_length
                subsequence   = sequence[start_position:stop_position]
            elif variant_type == 'frameshift_fusion':
                subsequence = sequence[start_position:]
                if subsequence.endswith('X'):
                    subsequence = subsequence[:-1]
            else:
                continue

            if '*' in subsequence:
                continue

            if 'X' in subsequence:
                continue

            if len(subsequence) < self.epitope_length:
                continue

            if subsequence in fasta_sequences:
                fasta_sequences[subsequence].append(line['index'])
            else:
                fasta_sequences[subsequence] = [line['index']]

        writer                  = open(self.output_file, 'w')
        key_writer              = open(self.output_key_file, 'w')
        count                   = 1
        for (subsequence, keys) in fasta_sequences.items():
            writer.writelines('>%s\n' % count)
            writer.writelines('%s\n' % subsequence)
            yaml.dump({count: keys}, key_writer, default_flow_style=False)
            count += 1

        reader.close()
        writer.close()
        key_writer.close()

