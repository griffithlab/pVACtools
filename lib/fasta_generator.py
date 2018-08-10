import csv
import re
import sys
from collections import OrderedDict, defaultdict
import yaml
from abc import ABCMeta
from Bio import SeqIO
import itertools
from lib.proximal_variant import ProximalVariant

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
        self.peptide_sequence_length    = kwargs['peptide_sequence_length']
        self.epitope_length             = kwargs['epitope_length']
        self.output_file                = kwargs['output_file']
        self.output_key_file            = kwargs['output_key_file']
        self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', None)
        self.proximal_variants_file     = kwargs.pop('proximal_variants_file', None)
        self.proximal_variants          = self.parse_proximal_variants_file()

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
        ##clip by wt sequence length, otherwise with deletions peptide_sequence_length may exceeds full wt sequence length,
        ##and the code below tries extracting ranges beyond the wt sequence
        peptide_sequence_length = min(2 * one_flanking_sequence_length + wildtype_amino_acid_length,len(full_wildtype_sequence))

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
        return start_position, wildtype_subsequence, mutation_start_subsequence

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

                proximal_variant_position = int(protein_position) - 1 - mutation_offset
                if proximal_variant_position <= 0 or proximal_variant_position >= len(wildtype_subsequence):
                    continue
                if len(proximal_variant_wildtype_amino_acid) != len(proximal_variant_mutant_amino_acid):
                    print("Nearby variant is not a missense mutation. Skipping.")
                    continue
                if wildtype_subsequence[proximal_variant_position] != proximal_variant_wildtype_amino_acid:
                    sys.exit(
                        "Error when processing proximal variant.\n" +
                        "The wildtype amino acid for variant %s with substring %s is different than expected.\n" % (somatic_variant_index, wildtype_subsequence) +
                        "Actual wildtype amino acid: %s\n" % wildtype_subsequence[proximal_variant_position] +
                        "Wildtype amino acid of the proximal_variant: %s" % proximal_variant_wildtype_amino_acid
                    )
                wildtype_subsequence_with_proximal_variants = wildtype_subsequence_with_proximal_variants[:proximal_variant_position] + proximal_variant_mutant_amino_acid + wildtype_subsequence_with_proximal_variants[proximal_variant_position+1:]
        return wildtype_subsequence_with_proximal_variants

    def execute(self):
        peptide_sequence_length = self.peptide_sequence_length
        reader                  = open(self.input_file, 'r')
        tsvin                   = csv.DictReader(reader, delimiter='\t')
        fasta_sequences         = OrderedDict()
        for line in tsvin:
            variant_type = line['variant_type']
            full_wildtype_sequence = line['wildtype_amino_acid_sequence']
            if variant_type == 'FS':
                if line['amino_acid_change'] is not None and line['amino_acid_change'].split('/')[0] == '-':
                    position = int(line['protein_position'].split('-', 1)[0])
                else:
                    position = int(line['protein_position'].split('-', 1)[0]) - 1
            elif variant_type == 'missense' or variant_type == 'inframe_ins':
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
                mutation_start_position, wildtype_subsequence, left_flanking_subsequence = self.get_frameshift_subsequences(position, full_wildtype_sequence, peptide_sequence_length, line)
                downstream_sequence = line['downstream_amino_acid_sequence']
                if self.downstream_sequence_length and len(downstream_sequence) > self.downstream_sequence_length:
                    downstream_sequence = downstream_sequence[0:self.downstream_sequence_length]
                wildtype_subsequence = self.add_proximal_variants(line['index'], wildtype_subsequence, mutation_start_position, position, True)
                left_flanking_subsequence_with_proximal_variants = self.add_proximal_variants(line['index'], left_flanking_subsequence, mutation_start_position, position, False)
                #The caveat here is that if a nearby variant is in the downstream sequence, the protein sequence would be further altered, which we aren't taking into account.
                #we would need to recalculate the downstream protein sequence taking all downstream variants into account.
                mutant_subsequence = left_flanking_subsequence_with_proximal_variants + downstream_sequence
            else:
                mutation_start_position, wildtype_subsequence = self.get_wildtype_subsequence(position, full_wildtype_sequence, wildtype_amino_acid_length, peptide_sequence_length, line)
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
        peptide_sequence_length = self.peptide_sequence_length
        reader                  = open(self.input_file, 'r')
        tsvin                   = csv.DictReader(reader, delimiter='\t')
        fasta_sequences         = OrderedDict()
        for line in tsvin:
            variant_type = line['variant_type']
            position     = int(line['protein_position'])
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
            else:
                continue

            if subsequence.endswith('X'):
                subsequence = subsequence[:-1]

            if '*' in subsequence:
                continue

            if 'X' in subsequence:
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
        self.input_file       = kwargs['input_file']
        self.output_file      = kwargs['output_file']
        self.output_key_file  = kwargs['output_key_file']

    def execute(self):
        seq_dict = dict()
        for record in SeqIO.parse(self.input_file, "fasta"):
            seq_dict[record.id] = str(record.seq)
        seq_keys = sorted(seq_dict)

        seq_tuples = list(itertools.permutations(seq_keys, 2))
        epitopes = dict()
        rev_lookup = dict()

        fasta_sequences = OrderedDict()
        for comb in seq_tuples:
            seq1 = comb[0]
            seq2 = comb[1]
            for length in range(8, 11):
                seq_ID = seq1 + "|" + seq2
                trunc_seq1 = seq_dict[seq1][(len(seq_dict[seq1]) - length):len(seq_dict[seq1])]
                trunc_seq2 = seq_dict[seq2][0:(length - 1)]
                epitopes[seq_ID] = trunc_seq1 + trunc_seq2
                rev_lookup[(trunc_seq1 + trunc_seq2)] = seq_ID

                spacers = ["HH", "HHC", "HHH", "HHHD", "HHHC", "AAY", "HHHH", "HHAA", "HHL", "AAL"]
                for this_spacer in spacers:
                    seq_ID = seq1 + "|" + this_spacer + "|" + seq2
                    epitopes[seq_ID] = (trunc_seq1 + this_spacer + trunc_seq2)
                    rev_lookup[(trunc_seq1 + this_spacer + trunc_seq2)] = seq_ID

        for seq_id in epitopes:
            sequence = epitopes[seq_id]
            fasta_sequences.setdefault(sequence, []).append(seq_id)

        writer = open(self.output_file, 'w')
        key_writer = open(self.output_key_file, 'w')
        count  = 1
        for (subsequence, keys) in sorted(fasta_sequences.items()):
            writer.writelines('>%s\n' % count)
            writer.writelines('%s\n' % subsequence)
            yaml.dump({count: keys}, key_writer, default_flow_style=False)
            count += 1

        writer.close()
        key_writer.close()
        self.epitopes = epitopes
        self.seq_tuples = seq_tuples
        self.seq_dict = seq_dict
        self.seq_keys = seq_keys
