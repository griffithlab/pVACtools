import os
import csv
import re
import tempfile
from .input_file_converter import *
from .fasta_generator import *
import math

class PvacvectorInputFastaGenerator():
    def __init__(self, pvacseq_tsv, input_vcf, output_dir, n_mer):
        self.input_tsv = pvacseq_tsv
        self.input_vcf = input_vcf
        self.output_dir = output_dir
        self.output_file = os.path.join(self.output_dir, "vector_input.fa")
        self.n_mer = int(n_mer)

    def parse_choosen_epitopes(self):
        epitopes = {}
        with open(self.input_tsv, 'r') as input_f:
            reader = csv.DictReader(input_f, delimiter = "\t")
            for line in reader:
                consequence = line['Variant Type']
                mutation = line['Mutation']
                pos = line['Protein Position']
                subpeptide_pos = line['Sub-peptide Position']
                gene_name = line['Gene Name']
                mt_epitope_seq = line['MT Epitope Seq']
                transcript = line['Transcript']

                if consequence == 'FS':
                    amino_acid_change_position = "%s%s/%s" % (pos, line['Reference'], line['Variant'])
                else:
                    amino_acid_change_position = pos + mutation
                index = '%s.%s.%s.%s_pos%s_len%s' % (gene_name, transcript, consequence, amino_acid_change_position, line['Sub-peptide Position'], line['Peptide Length'])
                epitopes[index] = mt_epitope_seq
        return epitopes

    #get necessary data from initial pvacseq input vcf
    def parse_original_vcf(self):
        tsv_file = tempfile.NamedTemporaryFile()
        VcfConverter(**{'input_file': self.input_vcf, 'output_file': tsv_file.name}).execute()
        fasta_file = tempfile.NamedTemporaryFile()
        key_file = tempfile.NamedTemporaryFile()
        FastaGenerator(**{
            'input_file': tsv_file.name,
            'peptide_sequence_length': self.n_mer + 2 * 8,
            'epitope_length': 8,
            'output_file': fasta_file.name,
            'output_key_file': key_file.name,
        }).execute()

        with open(key_file.name, 'r') as fasta_key_file:
            keys = yaml.load(fasta_key_file)

        dataframe = OrderedDict()
        with open(fasta_file.name, 'r') as fasta_file:
            for line in fasta_file:
                key      = line.rstrip().replace(">","")
                sequence = fasta_file.readline().rstrip()
                if key.startswith('WT'):
                    continue
                ids      = keys[int(key)]
                for id in ids:
                    (type, index) = id.split('.', 1)
                    dataframe[index] = sequence
        return dataframe

    def write_output_fasta(self, peptides):
        with open(self.output_file, 'w') as out_f:
            for index, peptide in peptides.items():
                out_f.write(">%s\n" % index)
                out_f.write("%s\n" % peptide)
                print("ID: " + index + ", sequence: " + peptide)
        out_f.close()
        print("FASTA file written")

    def execute(self):
        epitopes = self.parse_choosen_epitopes()
        transcripts_dict = self.parse_original_vcf()
        extracted_peptides = self.extract_peptide_sequences(transcripts_dict, epitopes)
        self.write_output_fasta(extracted_peptides)

    def extract_peptide_sequences(self, transcript_dict, epitopes):
        extracted_peptides = {}
        for index, epitope in sorted(epitopes.items()):
            p = re.compile('^(.+)_pos([0-9]+)_len([0-9]+)$')
            identifier = p.search(index).group(1)
            original_position = p.search(index).group(2)
            length = p.search(index).group(3)
            full_sequence = transcript_dict[identifier]
            #find all occurrences of the epitope in the full sequence
            occurrences = [n for n in range(len(full_sequence)) if full_sequence.find(epitope, n) == n]
            #epitope occurs once in the full sequence
            if len(occurrences) == 1:
                new_position = occurrences[0]
            #epitope occurs multiple times
            else:
                #find the occurence that is closest to the original positon and use that
                new_position = min(occurences, key=lambda x:abs(original_position-1-x))

            length = int(length)
            n_mer = int(self.n_mer)
            if len(full_sequence) <= n_mer:
                extracted_sequence = full_sequence
            else:
                wingspan = (n_mer - length) / 2
                start = math.ceil(new_position - wingspan)
                end = math.ceil(new_position + length + wingspan)
                if start < 0:
                    remainder = math.fabs(start)
                    end += remainder
                if end > len(full_sequence):
                    remainder = end - len(full_sequence)
                    start -= remainder
                extracted_sequence = full_sequence[start:end]
            exists = False
            for other_index, other_extracted_sequence in extracted_peptides.items():
                if other_extracted_sequence == extracted_sequence:
                    print("Sequence for item %s is identical to sequence for item %s. Skipping %s" % (index, other_index, index))
                    exists = True
            if not exists:
                extracted_peptides[index] = extracted_sequence
        return extracted_peptides

