import csv
from vaxrank.manufacturability import ManufacturabilityScores
from vaxrank.vaccine_peptide import *
from Bio import SeqIO

class PvacpeptideVaccinePeptide(VaccinePeptide):
    def __new__(cls, peptide):
        return VaccinePeptideBase.__new__(
            cls,
            mutant_protein_fragment="",
            mutant_epitope_predictions="",
            wildtype_epitope_predictions="",
            mutant_epitope_score="",
            wildtype_epitope_score="",
            num_mutant_epitopes_to_keep="",
            manufacturability_scores=ManufacturabilityScores.from_amino_acids(peptide)
        )

class CalculateManufacturability:
    def __init__(self, input_file, output_file, file_type='pVACseq'):
        self.input_file = input_file
        self.output_file = output_file
        self.file_type = file_type

    def manufacturability_headers(self):
        return [
            'cterm_7mer_gravy_score',
            'max_7mer_gravy_score',
            'difficult_n_terminal_residue',
            'c_terminal_cysteine',
            'c_terminal_proline',
            'cysteine_count',
            'n_terminal_asparagine',
            'asparagine_proline_bond_count'
        ]

    def append_manufacturability_metrics(self, line, peptide):
        metrics = peptide.manufacturability_scores
        line['cterm_7mer_gravy_score'] = metrics.cterm_7mer_gravy_score
        line['max_7mer_gravy_score'] = metrics.max_7mer_gravy_score
        line['difficult_n_terminal_residue'] = metrics.difficult_n_terminal_residue
        line['c_terminal_cysteine'] = metrics.c_terminal_cysteine
        line['c_terminal_proline'] = metrics.c_terminal_proline
        line['cysteine_count'] = metrics.cysteine_count
        line['n_terminal_asparagine'] = metrics.n_terminal_asparagine
        line['asparagine_proline_bond_count'] = metrics.asparagine_proline_bond_count
        return line

    def execute(self):
        if self.file_type == 'fasta':
            with open(self.output_file, 'w') as output_fh:
                writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=['id', 'peptide_sequence'] + self.manufacturability_headers(), extrasaction='ignore')
                writer.writeheader()
                for record in SeqIO.parse(self.input_file, "fasta"):
                    seq_num = record.id
                    peptide = str(record.seq)
                    line = {
                        'id': seq_num,
                        'peptide_sequence': peptide
                    }
                    peptide = PvacpeptideVaccinePeptide(peptide)
                    line = self.append_manufacturability_metrics(line, peptide)
                    writer.writerow(line)
        else:
            with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
                reader = csv.DictReader(input_fh, delimiter = "\t")
                writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=reader.fieldnames + self.manufacturability_headers(), extrasaction='ignore')
                writer.writeheader()
                for line in reader:
                    if self.file_type == 'pVACbind':
                        peptide = PvacpeptideVaccinePeptide(line['Epitope Seq'])
                    else:
                        peptide = PvacpeptideVaccinePeptide(line['MT Epitope Seq'])
                    line = self.append_manufacturability_metrics(line, peptide)
                    writer.writerow(line)
