import csv
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import shutil
import re

class CalculateReferenceProteomeSimilarity:
    def __init__(self, input_file, input_fasta, output_file, peptide_sequence_length, match_length=8, species='human', file_type='pVACseq'):
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

    def get_peptides(self):
        records = list(SeqIO.parse(self.input_fasta, "fasta"))
        if self.file_type == 'pVACbind':
            pass
        else:
            records_dict = {x.id.replace('MT.', ''): str(x.seq) for x in filter(lambda x: x.id.startswith('MT.'), records)}
        return records_dict

    def extract_n_mer(self, full_peptide, subpeptide_position, mutation_position, mt_length):
        flanking_sequence_length = self.match_length - 1
        mt_start = (subpeptide_position-1) + (mutation_position-1)
        start = mt_start - flanking_sequence_length
        if start < 0:
            start = 0
        end = mt_start + mt_length + flanking_sequence_length
        return full_peptide[start:end]

    def extract_n_mer_from_fs(self, full_peptide, epitope, peptide_sequence_length, subpeptide_position):
        if peptide_sequence_length%2 == 0:
            flanking_sequence_length = int(peptide_sequence_length/2)
        else:
            flanking_sequence_length = int((peptide_sequence_length-1)/2)
        start = subpeptide_position - 1 - flanking_sequence_length
        if start < 0:
            start = 0
        end = subpeptide_position + len(epitope) + flanking_sequence_length
        return full_peptide[start:end]

    def metric_headers(self):
        return ['Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 'Transcript', 'Hit ID', 'Hit Definition', 'Query Sequence', 'Match Sequence', 'Match Start', 'Match Stop']

    def execute(self):
        if self.species not in self.species_to_organism:
            print("Species {} not supported for Reference Proteome Similarity search. Skipping.".format(self.species))
            shutil.copy(self.input_file, self.output_file)
            return

        records_dict = self.get_peptides()

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
                if line['Variant Type'] == 'FS':
                    peptide = self.extract_n_mer_from_fs(records_dict[line['Index']], epitope, self.peptide_sequence_length, int(line['Sub-peptide Position']))
                else:
                    mt_amino_acids = line['Mutation'].split('/')[1]
                    if mt_amino_acids == '-':
                        mt_amino_acids = ''
                    peptide = self.extract_n_mer(records_dict[line['Index']], int(line['Sub-peptide Position']), int(line['Mutation Position']), len(mt_amino_acids))
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
