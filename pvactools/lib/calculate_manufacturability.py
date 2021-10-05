import csv
from vaxrank.manufacturability import ManufacturabilityScores
from Bio import SeqIO

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

    def append_manufacturability_metrics(self, line, metrics):
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
                writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=['id', 'peptide_sequence'] + self.manufacturability_headers(), extrasaction='ignore', restval='NA')
                writer.writeheader()
                for record in SeqIO.parse(self.input_file, "fasta"):
                    seq_num = record.id
                    sequence = str(record.seq)
                    line = {
                        'id': seq_num,
                        'peptide_sequence': sequence
                    }
                    if len(sequence) >= 7:
                        scores = ManufacturabilityScores.from_amino_acids(sequence)
                        line = self.append_manufacturability_metrics(line, scores)
                    writer.writerow(line)
        else:
            with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
                reader = csv.DictReader(input_fh, delimiter = "\t")
                writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames=reader.fieldnames + self.manufacturability_headers(), extrasaction='ignore', restval='NA')
                writer.writeheader()
                for line in reader:
                    if self.file_type == 'pVACbind' or self.file_type == 'pVACfuse':
                        sequence = line['Epitope Seq']
                    else:
                        sequence = line['MT Epitope Seq']
                    if len(sequence) >= 7:
                        scores = ManufacturabilityScores.from_amino_acids(sequence)
                        line = self.append_manufacturability_metrics(line, scores)
                    writer.writerow(line)
