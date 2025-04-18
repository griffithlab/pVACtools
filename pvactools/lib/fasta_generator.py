import csv
import json
import re
import sys
from collections import OrderedDict, defaultdict
import yaml
from abc import ABCMeta
from Bio import SeqIO
import itertools
import logging

from pvactools.lib.proximal_variant import ProximalVariant

csv.field_size_limit(sys.maxsize)

class FastaGenerator(metaclass=ABCMeta):
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

    def __init__(self, **kwargs):
        self.input_file                 = kwargs['input_file']
        self.flanking_sequence_length   = kwargs['flanking_sequence_length']
        self.epitope_length             = kwargs['epitope_length']
        self.output_file                = kwargs['output_file']
        self.output_key_file            = kwargs['output_key_file']
        self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', None)
        self.proximal_variants_file     = kwargs.pop('proximal_variants_file', None)
        self.trim_invalid_characters    = kwargs.pop('trim_invalid_characters', False)
        self.proximal_variants          = self.parse_proximal_variants_file()

    def invalid_characters(self):
        return ['*', 'X', '?']

    def contains_invalid_characters(self, sequence):
        for character in self.invalid_characters():
            if character in sequence:
                return True
        return False

    def trim_sequence(self, sequence):
        for character in self.invalid_characters():
            while character in sequence:
                invalid_character_pos = sequence.index(character)
                if invalid_character_pos < (len(sequence) - invalid_character_pos):
                    sequence = sequence[invalid_character_pos+1:]
                else:
                    sequence = sequence[0:invalid_character_pos]
        return sequence

    def position_out_of_bounds(self, position, sequence):
        return position > len(sequence)-1

    #This subroutine is a bit funky but it was designed that way to mirror
    #distance_from_end to increase code readability from the caller's perspective
    def distance_from_start(self, position, string):
        return position

    def distance_from_end(self, position, string):
        return len(string) - 1 - position

    def get_wildtype_subsequence(self, position, full_wildtype_sequence, wildtype_amino_acid_length, line):
        ##clip by wt sequence length, otherwise with deletions peptide_sequence_length may exceeds full wt sequence length,
        ##and the code below tries extracting ranges beyond the wt sequence
        peptide_sequence_length = min(2 * self.flanking_sequence_length + wildtype_amino_acid_length, len(full_wildtype_sequence))

        # We want to extract a subset from full_wildtype_sequence that is
        # peptide_sequence_length long so that the position ends
        # up in the middle of the extracted sequence.
        # If the position is too far toward the beginning or end of
        # full_wildtype_sequence there aren't enough amino acids on one side
        # to achieve this.
        if self.distance_from_start(position, full_wildtype_sequence) < self.flanking_sequence_length:
            wildtype_subsequence = full_wildtype_sequence[:peptide_sequence_length]
            mutation_position = position
        elif self.distance_from_end(position, full_wildtype_sequence) < self.flanking_sequence_length:
            start_position = len(full_wildtype_sequence) - peptide_sequence_length
            wildtype_subsequence = full_wildtype_sequence[start_position:]
            mutation_position = peptide_sequence_length - self.distance_from_end(position, full_wildtype_sequence) - 1
        elif self.distance_from_start(position, full_wildtype_sequence) >= self.flanking_sequence_length and self.distance_from_end(position, full_wildtype_sequence) >= self.flanking_sequence_length:
            start_position = position - self.flanking_sequence_length
            end_position   = start_position + peptide_sequence_length
            wildtype_subsequence = full_wildtype_sequence[start_position:end_position]
            mutation_position = self.flanking_sequence_length
        else:
            sys.exit("ERROR: Something went wrong during the retrieval of the wildtype sequence at position(%s, %s, %s)" % line['chromsome_name'], line['start'], line['stop'])

        return mutation_position, wildtype_subsequence

    def get_frameshift_subsequences(self, position, full_wildtype_sequence, full_mutant_sequence):
        if position < self.flanking_sequence_length:
            start_position = 0
        else:
            start_position = position - self.flanking_sequence_length
        wildtype_subsequence_stop_position = position + self.flanking_sequence_length
        wildtype_subsequence = full_wildtype_sequence[start_position:wildtype_subsequence_stop_position]
        if self.downstream_sequence_length:
            mutant_subsequence_stop_position = position + self.downstream_sequence_length
            mutant_subsequence = full_mutant_sequence[start_position:mutant_subsequence_stop_position]
        else:
            mutant_subsequence = full_mutant_sequence[start_position:]
        left_flanking_sequence = full_mutant_sequence[start_position:position]
        return wildtype_subsequence, mutant_subsequence, left_flanking_sequence

    def add_proximal_variants(self, somatic_variant_index, wildtype_subsequence, mutation_position, original_position, germline_variants_only):
        mutation_offset = original_position - mutation_position
        wildtype_subsequence_with_proximal_variants = wildtype_subsequence
        if somatic_variant_index in self.proximal_variants.keys():
            for (protein_position, lines) in self.proximal_variants[somatic_variant_index].items():
                if protein_position == original_position:
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
                    proximal_variant_start_position = int(protein_position.split('-')[0]) - 1 - mutation_offset
                else:
                    proximal_variant_start_position = int(protein_position) - 1 - mutation_offset
                proximal_variant_end_position = proximal_variant_start_position + len(proximal_variant_mutant_amino_acid)

                if proximal_variant_end_position <= 0 or proximal_variant_start_position >= len(wildtype_subsequence):
                    continue
                if len(proximal_variant_wildtype_amino_acid) != len(proximal_variant_mutant_amino_acid):
                    print("Nearby variant is not a missense mutation. Skipping.")
                    continue

                if proximal_variant_end_position > len(wildtype_subsequence):
                    #the DNP extends past the end of the wildtype_subsequence
                    missing_amino_acids_count = proximal_variant_end_position - len(wildtype_subsequence)
                    #remove proximal variant amino acids after wildtype subsquence end
                    proximal_variant_wildtype_amino_acid = proximal_variant_wildtype_amino_acid[:len(proximal_variant_wildtype_amino_acid) - missing_amino_acids_count]
                    proximal_variant_mutant_amino_acid = proximal_variant_mutant_amino_acid[:len(proximal_variant_mutant_amino_acid) - missing_amino_acids_count]

                if proximal_variant_start_position < 0:
                    #DNP starts before the beginning of the wildtype subsequence
                    missing_amino_acids_count = abs(proximal_variant_start_position)
                    proximal_variant_start_position = 0
                    #remove proximal variant amino acids before wildtype subsequence start
                    proximal_variant_wildtype_amino_acid = proximal_variant_wildtype_amino_acid[missing_amino_acids_count:]
                    proximal_variant_mutant_amino_acid = proximal_variant_mutant_amino_acid[missing_amino_acids_count:]

                if wildtype_subsequence[proximal_variant_start_position:proximal_variant_end_position] != proximal_variant_wildtype_amino_acid:
                    sys.exit(
                        "Error when processing proximal variant.\n" +
                        "The wildtype amino acid for variant %s with substring %s is different than expected.\n" % (somatic_variant_index, wildtype_subsequence) +
                        "Actual wildtype amino acid: %s\n" % wildtype_subsequence[proximal_variant_start_position:proximal_variant_end_position] +
                        "Wildtype amino acid of the proximal_variant: %s" % proximal_variant_wildtype_amino_acid
                    )

                wildtype_subsequence_with_proximal_variants = wildtype_subsequence_with_proximal_variants[:proximal_variant_start_position] + proximal_variant_mutant_amino_acid + wildtype_subsequence_with_proximal_variants[proximal_variant_end_position:]
        return wildtype_subsequence_with_proximal_variants

    def execute(self):
        reader                  = open(self.input_file, 'r')
        tsvin                   = csv.DictReader(reader, delimiter='\t')
        fasta_sequences         = OrderedDict()
        for line in tsvin:
            variant_type = line['variant_type']
            full_wildtype_sequence = line['wildtype_amino_acid_sequence']
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

            if self.position_out_of_bounds(position, full_wildtype_sequence):
                continue

            if variant_type == 'missense' and line['index'] in self.proximal_variants and line['protein_position'] in self.proximal_variants[line['index']]:
                codon_changes = [ item['codon_change'] for item in self.proximal_variants[line['index']][line['protein_position']] ]
                codon_changes.append(line['codon_change'])
                mutant_amino_acid_with_proximal_variants = ProximalVariant.combine_conflicting_variants(codon_changes)
            elif variant_type != 'FS':
                mutant_amino_acid_with_proximal_variants = mutant_amino_acid

            if variant_type == 'FS':
                full_mutant_sequence = line['frameshift_amino_acid_sequence']
                wildtype_subsequence, mutant_subsequence, left_flanking_subsequence = self.get_frameshift_subsequences(position, full_wildtype_sequence, full_mutant_sequence)
                mutation_start_position = len(left_flanking_subsequence)
                wildtype_subsequence = self.add_proximal_variants(line['index'], wildtype_subsequence, mutation_start_position, position, True)
                left_flanking_subsequence_with_proximal_variants = self.add_proximal_variants(line['index'], left_flanking_subsequence, mutation_start_position, position, False)
                #The caveat here is that if a nearby variant is in the downstream sequence, the protein sequence would be further altered, which we aren't taking into account.
                #we would need to recalculate the downstream protein sequence taking all downstream variants into account.
                mutant_subsequence = re.sub('^%s' % left_flanking_subsequence, left_flanking_subsequence_with_proximal_variants, mutant_subsequence)
            else:
                if variant_type == 'inframe_ins':
                    mutation_start_position, wildtype_subsequence = self.get_wildtype_subsequence(position, full_wildtype_sequence, len(mutant_amino_acid), line)
                else:
                    mutation_start_position, wildtype_subsequence = self.get_wildtype_subsequence(position, full_wildtype_sequence, wildtype_amino_acid_length, line)
                mutation_end_position = mutation_start_position + wildtype_amino_acid_length
                if wildtype_amino_acid != '-' and wildtype_amino_acid != wildtype_subsequence[mutation_start_position:mutation_end_position]:
                    if line['amino_acid_change'].split('/')[0].count('*') > 1:
                        print("Warning: Amino acid change is not sane - contains multiple stops. Skipping entry {}".format(line['index']))
                        continue
                    else:
                        sys.exit("ERROR: There was a mismatch between the actual wildtype amino acid sequence ({}) and the expected amino acid sequence ({}). Did you use the same reference build version for VEP that you used for creating the VCF?\n{}".format(wildtype_subsequence[mutation_start_position:mutation_end_position], wildtype_amino_acid, line))
                wildtype_subsequence_with_proximal_variants = self.add_proximal_variants(line['index'], wildtype_subsequence, mutation_start_position, position, False)
                wildtype_subsequence = self.add_proximal_variants(line['index'], wildtype_subsequence, mutation_start_position, position, True)
                if stop_codon_added:
                    mutant_subsequence = wildtype_subsequence_with_proximal_variants[:mutation_start_position] + mutant_amino_acid_with_proximal_variants
                else:
                    mutant_subsequence = wildtype_subsequence_with_proximal_variants[:mutation_start_position] + mutant_amino_acid_with_proximal_variants + wildtype_subsequence_with_proximal_variants[mutation_end_position:]

            if '*' in wildtype_subsequence or '*' in mutant_subsequence:
                continue

            if 'X' in wildtype_subsequence or 'X' in mutant_subsequence:
                continue

            if 'U' in wildtype_subsequence or 'U' in mutant_subsequence:
                print("Warning. Sequence contains unsupported amino acid U. Skipping entry {}".format(line['index']))
                continue

            if mutant_subsequence in wildtype_subsequence:
                #This is not a novel peptide
                continue

            if len(wildtype_subsequence) < self.epitope_length or len(mutant_subsequence) < self.epitope_length:
                continue

            variant_id = line['index']
            for designation, subsequence in zip(['WT', 'MT'], [wildtype_subsequence, mutant_subsequence]):
                key = '%s.%s' % (designation, variant_id)
                fasta_sequences.setdefault(subsequence, []).append(key)

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
        reader                  = open(self.input_file, 'r')
        tsvin                   = csv.DictReader(reader, delimiter='\t')
        fasta_sequences         = OrderedDict()
        for line in tsvin:
            variant_type = line['variant_type']
            position     = int(line['protein_position'])
            sequence     = line['fusion_amino_acid_sequence']
            if position < self.flanking_sequence_length:
                start_position = 0
            else:
                start_position = position - self.flanking_sequence_length

            if variant_type == 'inframe_fusion':
                stop_position = position + self.flanking_sequence_length
                subsequence   = sequence[start_position:stop_position]
            elif variant_type == 'frameshift_fusion':
                if self.downstream_sequence_length:
                    stop_position = position + self.downstream_sequence_length
                    subsequence = sequence[start_position:stop_position]
                else:
                    subsequence = sequence[start_position:]
            else:
                continue

            if subsequence.endswith('X'):
                subsequence = subsequence[:-1]

            if self.contains_invalid_characters(subsequence):
                if self.trim_invalid_characters:
                    subsequence = self.trim_sequence(subsequence)
                else:
                    continue

            if len(subsequence) < self.epitope_length:
                continue

            fasta_sequences.setdefault(subsequence, []).append(line['index'])

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

class VectorFastaGenerator():
    def __init__(self, **kwargs):
        self.input_file         = kwargs['input_file']
        self.output_file_prefix = kwargs['output_file_prefix']
        self.epitope_lengths    = kwargs['epitope_lengths']
        self.spacer             = kwargs['spacer']
        self.junctions_to_test  = kwargs['junctions_to_test']
        self.clip_length        = kwargs['clip_length']
        self.output_files = []

    def execute(self):
        seq_dict = dict()
        best_peptide = dict()
        for record in SeqIO.parse(self.input_file, "fasta"):
            seq_dict[record.id] = str(record.seq)
            description = record.description.replace("{} ".format(record.id), "")
            if description != "":
               try:
                   best_peptide[record.id] = json.loads(description)['Best Peptide']
               except:
                   pass

        for length in self.epitope_lengths:
            epitopes = dict()
            fasta_sequences = OrderedDict()
            wingspan_length = length - 1
            warnings = set()
            for (seq1, seq2) in self.junctions_to_test:
                seq1_seq = seq_dict[seq1]
                seq2_seq = seq_dict[seq2]
                for left_clip_length in range(0, self.clip_length+1):
                    for right_clip_length in range(0, self.clip_length+1):
                        #These combinations would've already been tested in previous attempts with lower clip lengths and can be skipped
                        if left_clip_length < self.clip_length and right_clip_length < self.clip_length:
                            continue
                        if seq1 in best_peptide:
                            seq1_best_peptide = best_peptide[seq1]
                            last_position = seq1_seq.rindex(seq1_best_peptide) + len(seq1_best_peptide)
                            end_distance = len(seq1_seq) - last_position
                            if left_clip_length > end_distance:
                                warnings.add("Clipping {} amino acids off the end of peptide {} would clip the best peptide. Skipping.".format(left_clip_length, seq1))
                                continue
                        if seq2 in best_peptide:
                            seq2_best_peptide = best_peptide[seq2]
                            first_position = seq2_seq.index(seq2_best_peptide)
                            if right_clip_length > first_position:
                                warnings.add("Clipping {} amino acids off the start of peptide {} would clip the best peptide. Skipping.".format(right_clip_length, seq2))
                                continue
                        trunc_seq1 = seq1_seq[(len(seq1_seq) - wingspan_length - left_clip_length):(len(seq1_seq) - left_clip_length)]
                        trunc_seq2 = seq2_seq[(0 + right_clip_length):wingspan_length + right_clip_length]

                        if self.spacer != 'None':
                            seq_ID = "{}|{}|{}|{}|{}".format(seq1, left_clip_length, self.spacer, right_clip_length, seq2)
                            epitopes[seq_ID] = (trunc_seq1 + self.spacer + trunc_seq2)
                        else:
                            seq_ID = "{}|{}|{}|{}".format(seq1, left_clip_length, right_clip_length, seq2)
                            epitopes[seq_ID] = trunc_seq1 + trunc_seq2
            for warning in list(warnings):
                logging.info(warning)

            for seq_id in epitopes:
                sequence = epitopes[seq_id]
                if len(sequence) < length:
                    continue
                fasta_sequences.setdefault(sequence, []).append(seq_id)

            output_file = "{}.{}.tsv".format(self.output_file_prefix, length)
            self.output_files.append(output_file)
            output_key_file = "{}.key".format(output_file)
            writer = open(output_file, 'w')
            key_writer = open(output_key_file, 'w')
            count  = 1
            for (subsequence, keys) in sorted(fasta_sequences.items()):
                writer.writelines('>%s\n' % count)
                writer.writelines('%s\n' % subsequence)
                yaml.dump({count: keys}, key_writer, default_flow_style=False)
                count += 1

            writer.close()
            key_writer.close()
