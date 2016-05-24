import argparse
import csv
import re
import sys

#This subroutine is a bit funky but it was designed that way to mirror
#distance_from_end to increase code readability from the caller's perspective
def distance_from_start(position, string):
    return position

def distance_from_end(position, string):
    return len(string) - 1 - position;

def get_wildtype_subsequence_for_printing(position, wildtype_sequence, peptide_sequence_length, line):
    chromosome = line[0]
    start      = line[1]
    stop       = line[2]
    #If the wildtype sequence is shorter than the desired peptide sequence
    #length we use the wildtype sequence length instead so that the extraction
    #algorithm below works correctly
    if len(wildtype_sequence) < peptide_sequence_length:
        peptide_sequence_length = len(wildtype_sequence)
        print("Wildtype sequence length is shorter than desired peptide sequence length at position (%s, %s, %s). Using wildtype sequence length (%s) instead." % chromosome, start, stop, peptide_sequence_length)

    # We want to extract a subset from @arr_wildtype_sequence that is
    # $peptide_sequence_length long so that the $position ends
    # up in the middle of the extracted sequence.
    # If the $position is too far toward the beginning or end of
    # @arr_wildtype_sequence there aren't enough amino acids on one side
    # to achieve this.
    one_flanking_sequence_length = int((peptide_sequence_length - 1) / 2)
    if distance_from_start(position, wildtype_sequence) < one_flanking_sequence_length:
        wildtype_subsequence = wildtype_sequence[0 : peptide_sequence_length]
        mutation_position = position
    elif distance_from_end(position, wildtype_sequence) < one_flanking_sequence_length:
        wildtype_subsequence = wildtype_sequence[len(wildtype_sequence)-peptide_sequence_length : len(wildtype_sequence)]
        mutation_position = peptide_sequence_length - distance_from_end(position, wildtype_sequence) - 1
    elif distance_from_start(position, wildtype_sequence) >= one_flanking_sequence_length and distance_from_end(position, wildtype_sequence) >= one_flanking_sequence_length:
        wildtype_subsequence = wildtype_sequence[position-one_flanking_sequence_length : position+one_flanking_sequence_length+1]
        mutation_position = one_flanking_sequence_length
    else:
        sys.exit("Something went wrong during the retrieval of the wildtype sequence at position(%s, %s, %s)" % chromsome, start, stop)

    return mutation_position, wildtype_subsequence

parser = argparse.ArgumentParser(description='')
parser.add_argument('input_file', type=argparse.FileType('r'), help='input list of variants',)
parser.add_argument('peptide_sequence_length', type=int, help='length of the peptide sequence')
parser.add_argument('output_file', type=argparse.FileType('w'), help='output FASTA file')

args = parser.parse_args()

tsvin = csv.reader(args.input_file, delimiter='\t')

for line in tsvin:
    pattern = re.compile('([A-Z])(\d+)([A-Z])');
    match = pattern.match(line[7]);
    if match is not None:
        wildtype_amino_acid, position, mutant_amino_acid = match.group(1, 2, 3)
        position = int(position) - 1

        wildtype_sequence = line[9];
        if wildtype_amino_acid != wildtype_sequence[position]:
            continue
        else:
            mutation_position, wildtype_subsequence = get_wildtype_subsequence_for_printing(position, wildtype_sequence, args.peptide_sequence_length, line)
            mutant_subsequence = list(wildtype_subsequence)

            mutant_subsequence[mutation_position] = mutant_amino_acid
            if len(wildtype_subsequence) > 0:
                for designation, subsequence in zip(['WT', 'MT'], [wildtype_subsequence, ''.join(mutant_subsequence)]):
                    fasta_header = '.'.join(['>' + designation, line[5], line[7]])
                    args.output_file.writelines([fasta_header + "\n", subsequence + "\n"])
            else:
                print(join("\t", 'NULL', position))

args.input_file.close()
args.output_file.close()

