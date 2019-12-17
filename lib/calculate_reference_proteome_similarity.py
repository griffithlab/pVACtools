import csv
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import shutil
import re
import os

class CalculateReferenceProteomeSimilarity:
    def __init__(self, input_file, input_fasta, output_file, peptide_sequence_length, match_length=8, species='human', file_type='vcf'):
        self.input_file = input_file
        self.input_fasta = input_fasta
        self.output_file = output_file
        self.metric_file = "{}.reference_matches".format(output_file)
        self.peptide_sequence_length = peptide_sequence_length
        self.match_length = match_length
        self.species = species
        self.file_type = file_type
        self.species_to_organism = {
            'human': 'Homo sapiens',
            'mouse': 'Mus musculus',
            'chimpanzee': 'Pan troglodytes',
            'macaque': 'Macaca',
            'cow': 'Bos taurus',
            'horse': 'Equus caballus',
            'pig': 'Sus scrofa',
        }

    def reference_match_headers(self):
        return [
            'Reference Match',
        ]

    def get_mt_peptides(self):
        records = list(SeqIO.parse(self.input_fasta, "fasta"))
        if self.file_type == 'vcf':
            records_dict = {x.id.replace('MT.', ''): str(x.seq) for x in filter(lambda x: x.id.startswith('MT.'), records)}
        elif self.file_type == 'bedpe':
            records_dict = {x.id: str(x.seq) for x in records}
        return records_dict

    def get_wt_peptides(self):
        if self.file_type == 'vcf':
            records = list(SeqIO.parse(self.input_fasta, "fasta"))
            records_dict = {x.id.replace('WT.', ''): str(x.seq) for x in filter(lambda x: x.id.startswith('WT.'), records)}
        else:
            return {}
        return records_dict

    def extract_n_mer(self, full_peptide, subpeptide_position, mutation_position, mt_length):
        #For non-frameshifts this ensures that we only test match_length epitopes that overlap the mutation
        #If we extract a larger region, we will get false-positive matches against the reference proteome
        #from the native wildtype portion of the peptide
        flanking_sequence_length = self.match_length - 1
        mt_start = (subpeptide_position-1) + (mutation_position-1)
        start = mt_start - flanking_sequence_length
        if start < 0:
            start = 0
        end = mt_start + mt_length + flanking_sequence_length
        return full_peptide[start:end]

    def extract_n_mer_from_fs(self, full_peptide, wt_peptide, epitope, peptide_sequence_length, subpeptide_position):
        #For frameshifts we want to test all downstream epitopes that would be part of the peptide_sequence_length
        #peptide since they are all potentially novel
        if peptide_sequence_length%2 == 0:
            flanking_sequence_length = int(peptide_sequence_length/2)
        else:
            flanking_sequence_length = int((peptide_sequence_length-1)/2)
        start = subpeptide_position - 1 - flanking_sequence_length
        if start < 0:
            start = 0
        #This catches cases where the start position would cause too many leading wildtype amino acids, which would result
        #in false-positive reference matches
        diff_position = [i for i in range(len(wt_peptide)) if wt_peptide[i] != full_peptide[i]][0]
        min_start = diff_position - self.match_length + 1 
        if min_start > start:
            start = min_start
        end = start + flanking_sequence_length + len(epitope) + flanking_sequence_length
        return full_peptide[start:end]

    def metric_headers(self):
        return ['Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 'Transcript', 'Peptide', 'Hit ID', 'Hit Definition', 'Query Sequence', 'Match Sequence', 'Match Start', 'Match Stop']

    def execute(self):
        if self.species not in self.species_to_organism:
            print("Species {} not supported for Reference Proteome Similarity search. Skipping.".format(self.species))
            shutil.copy(self.input_file, self.output_file)
            return

        mt_records_dict = self.get_mt_peptides()
        wt_records_dict = self.get_wt_peptides()

        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh, open(self.metric_file, 'w') as metric_fh:
            reader = csv.DictReader(input_fh, delimiter="\t")
            writer = csv.DictWriter(output_fh, delimiter="\t", fieldnames=reader.fieldnames + self.reference_match_headers(), extrasaction='ignore')
            metric_writer = csv.DictWriter(metric_fh, delimiter="\t", fieldnames=self.metric_headers(), extrasaction='ignore')
            writer.writeheader()
            metric_writer.writeheader()
            for line in reader:
                if self.file_type == 'pVACbind':
                    epitope = line['Epitope Seq']
                else:
                    epitope = line['MT Epitope Seq']
                if self.file_type == 'vcf':
                    if line['Variant Type'] == 'FS':
                        peptide = self.extract_n_mer_from_fs(mt_records_dict[line['Index']], wt_records_dict[line['Index']], epitope, self.peptide_sequence_length, int(line['Sub-peptide Position']))
                    else:
                        mt_amino_acids = line['Mutation'].split('/')[1]
                        if mt_amino_acids == '-':
                            mt_amino_acids = ''
                        peptide = self.extract_n_mer(mt_records_dict[line['Index']], int(line['Sub-peptide Position']), int(line['Mutation Position']), len(mt_amino_acids))
                elif self.file_type == 'bedpe':
                    peptide = mt_records_dict[line['Index']]
                else:
                result_handle = NCBIWWW.qblast("blastp", "refseq_protein", peptide, entrez_query="{} [Organism]".format(self.species_to_organism[self.species]))
                reference_match = False
                for blast_record in NCBIXML.parse(result_handle):
                    if len(blast_record.alignments) > 0:
                        for alignment in blast_record.alignments:
                            for hsp in alignment.hsps:
                                matches = re.split('\+| ', hsp.match)
                                for match in matches:
                                    if len(match) >= self.match_length:
                                        reference_match = True
                                        metric_line = line.copy()
                                        metric_line['Peptide'] = peptide
                                        metric_line['Hit ID'] = alignment.hit_id
                                        metric_line['Hit Definition'] = alignment.hit_def
                                        metric_line['Query Sequence'] = hsp.query
                                        metric_line['Match Sequence'] = hsp.match
                                        metric_line['Match Start'] = hsp.sbjct_start
                                        metric_line['Match Stop'] = hsp.sbjct_end
                                        metric_writer.writerow(metric_line)
                                        break
                line['Reference Match'] = reference_match
                writer.writerow(line)
