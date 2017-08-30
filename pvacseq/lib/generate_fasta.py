import argparse
import csv
import re
import sys
import yaml
from collections import OrderedDict

csv.field_size_limit(sys.maxsize)

def position_out_of_bounds(position, sequence):
    return position > len(sequence)-1

#This subroutine is a bit funky but it was designed that way to mirror
#distance_from_end to increase code readability from the caller's perspective
def distance_from_start(position, string):
    return position

def distance_from_end(position, string):
    return len(string) - 1 - position

def determine_peptide_sequence_length(full_wildtype_sequence_length, peptide_sequence_length, line):
    actual_peptide_sequence_length = peptide_sequence_length

    #If the wildtype sequence is shorter than the desired peptide sequence
    #length we use the wildtype sequence length instead so that the extraction
    #algorithm below works correctly
    if full_wildtype_sequence_length < actual_peptide_sequence_length:
        actual_peptide_sequence_length = full_wildtype_sequence_length
        print('Wildtype sequence length is shorter than desired peptide sequence length at position (%s, %s, %s). Using wildtype sequence length (%s) instead.' % (line['chromosome_name'], line['start'], line['stop'], actual_peptide_sequence_length))

    return actual_peptide_sequence_length

def determine_flanking_sequence_length(full_wildtype_sequence_length, peptide_sequence_length, line):
    actual_peptide_sequence_length = determine_peptide_sequence_length(full_wildtype_sequence_length, peptide_sequence_length, line)
    if actual_peptide_sequence_length%2 == 0:
        return (actual_peptide_sequence_length-2) / 2
    else:
        return (actual_peptide_sequence_length-1) / 2

def get_wildtype_subsequence(position, full_wildtype_sequence, wildtype_amino_acid_length, peptide_sequence_length, line):
    one_flanking_sequence_length = int(determine_flanking_sequence_length(len(full_wildtype_sequence), peptide_sequence_length, line))
    peptide_sequence_length = 2 * one_flanking_sequence_length + wildtype_amino_acid_length

    # We want to extract a subset from full_wildtype_sequence that is
    # peptide_sequence_length long so that the position ends
    # up in the middle of the extracted sequence.
    # If the position is too far toward the beginning or end of
    # full_wildtype_sequence there aren't enough amino acids on one side
    # to achieve this.
    if distance_from_start(position, full_wildtype_sequence) < one_flanking_sequence_length:
        wildtype_subsequence = full_wildtype_sequence[:peptide_sequence_length]
        mutation_position = position
    elif distance_from_end(position, full_wildtype_sequence) < one_flanking_sequence_length:
        start_position = len(full_wildtype_sequence) - peptide_sequence_length
        wildtype_subsequence = full_wildtype_sequence[start_position:]
        mutation_position = peptide_sequence_length - distance_from_end(position, full_wildtype_sequence) - 1
    elif distance_from_start(position, full_wildtype_sequence) >= one_flanking_sequence_length and distance_from_end(position, full_wildtype_sequence) >= one_flanking_sequence_length:
        start_position = position - one_flanking_sequence_length
        end_position   = start_position + peptide_sequence_length
        wildtype_subsequence = full_wildtype_sequence[start_position:end_position]
        mutation_position = one_flanking_sequence_length
    else:
        sys.exit("ERROR: Something went wrong during the retrieval of the wildtype sequence at position(%s, %s, %s)" % line['chromsome_name'], line['start'], line['stop'])

    return mutation_position, wildtype_subsequence

def get_frameshift_subsequences(position, full_wildtype_sequence, peptide_sequence_length, line):
    one_flanking_sequence_length = determine_flanking_sequence_length(len(full_wildtype_sequence), peptide_sequence_length, line)
    if position < one_flanking_sequence_length:
        start_position = 0
    else:
        start_position = int(position - one_flanking_sequence_length)
    wildtype_subsequence_stop_position = int(position + one_flanking_sequence_length)
    mutation_subsequence_stop_position = int(position)
    wildtype_subsequence = full_wildtype_sequence[start_position:wildtype_subsequence_stop_position]
    mutation_start_subsequence = full_wildtype_sequence[start_position:mutation_subsequence_stop_position]
    return wildtype_subsequence, mutation_start_subsequence

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq generate_fasta')
    parser.add_argument('input_file', type=argparse.FileType('r'), help='input list of variants',)
    parser.add_argument('peptide_sequence_length', type=int, help='length of the peptide sequence')
    parser.add_argument('epitope_length', type=int, help='length of subpeptides(epitopes) to predict')
    parser.add_argument('output_file', type=argparse.FileType('w'), help='output FASTA file')
    parser.add_argument('output_key_file', type=argparse.FileType('w'), help='output FASTA key file')
    parser.add_argument("-d", "--downstream-sequence-length", type=int, help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file.")
    args = parser.parse_args(args_input)

    peptide_sequence_length = args.peptide_sequence_length
    tsvin                   = csv.DictReader(args.input_file, delimiter='\t')
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

        if position_out_of_bounds(position, full_wildtype_sequence):
            continue

        if variant_type == 'FS':
            wildtype_subsequence, mutant_subsequence = get_frameshift_subsequences(position, full_wildtype_sequence, peptide_sequence_length, line)
            downstream_sequence = line['downstream_amino_acid_sequence']

            if args.downstream_sequence_length and len(downstream_sequence) > args.downstream_sequence_length:
                downstream_sequence = downstream_sequence[0:args.downstream_sequence_length]
            mutant_subsequence += downstream_sequence
        else:
            mutation_start_position, wildtype_subsequence = get_wildtype_subsequence(position, full_wildtype_sequence, wildtype_amino_acid_length, peptide_sequence_length, line)
            mutation_end_position = mutation_start_position + wildtype_amino_acid_length
            if wildtype_amino_acid != '-' and wildtype_amino_acid != wildtype_subsequence[mutation_start_position:mutation_end_position]:
                sys.exit("ERROR: There was a mismatch between the actual wildtype amino acid and the expected amino acid. Did you use the same reference build version for VEP that you used for creating the VCF?\n%s" % line)
            mutant_subsequence = wildtype_subsequence[:mutation_start_position] + mutant_amino_acid + wildtype_subsequence[mutation_end_position:]

        if '*' in wildtype_subsequence or '*' in mutant_subsequence:
            continue

        if 'X' in wildtype_subsequence or 'X' in mutant_subsequence:
            continue

        if len(wildtype_subsequence) < args.epitope_length or len(mutant_subsequence) < args.epitope_length:
            continue

        variant_id = line['index']
        for designation, subsequence in zip(['WT', 'MT'], [wildtype_subsequence, mutant_subsequence]):
            key = '%s.%s' % (designation, variant_id)
            if subsequence in fasta_sequences:
                fasta_sequences[subsequence].append(key)
            else:
                fasta_sequences[subsequence] = [key]

    count = 1
    for (subsequence, keys) in fasta_sequences.items():
        args.output_file.writelines('>%s\n' % count)
        args.output_file.writelines('%s\n' % subsequence)
        yaml.dump({count: keys}, args.output_key_file, default_flow_style=False)
        count += 1

    args.input_file.close()
    args.output_file.close()
    args.output_key_file.close()

if __name__ == '__main__':
    main()
