import csv
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import shutil
import re
import os
from collections import defaultdict
from time import sleep

class CalculateReferenceProteomeSimilarity:
    def __init__(self, input_file, input_fasta, output_file, match_length=8, species='human', file_type='vcf'):
        self.input_file = input_file
        self.input_fasta = input_fasta
        self.output_file = output_file
        self.metric_file = "{}.reference_matches".format(output_file)
        self.match_length = match_length
        self.species = species
        self.file_type = file_type
        self.species_to_organism = {
            'human': 'Homo sapiens',
            'atlantic salmon': 'Salmo salar',
            'black-headed spider monkey': 'Ateles fusciceps',
            'blue monkey': 'Cercopithecus mitis',
            'bonobo': 'Pan paniscus',
            'bornean orangutan': 'Pongo pygmaeus',
            'brown-mantled tamarin': 'Saguinus fuscicollis',
            'chimpanzee': 'Pan troglodytes',
            'common marmoset': 'Callithrix jacchus',
            'common squirrel monkey': 'Saimiri sciureus',
            'cottontop tamarin': 'Saguinus oedipus',
            'cow': 'Bos taurus',
            'crab-eating macaque': 'Macaca fascicularis',
            'dog': 'Canis lupus familiaris',
            "Geoffroy's tamarin": 'Saguinus geoffroyi',
            'golden lion tamarin': 'Leontopithecus rosalia',
            'gorilla': 'Gorilla gorilla',
            'grivet': 'Chlorocebus aethiops',
            'hamadryas baboon': 'Papio hamadryas',
            'horse': 'Equus caballus',
            'lar gibbon': 'Hylobates lar',
            'mouse': 'Mus musculus',
            'moustached tamarin': 'Saguinus mystax',
            'olive baboon': 'Papio anubis',
            'pig': 'Sus scrofa',
            'rainbow trout': 'Oncorhynchus mykiss',
            'rhesus macaque': 'Macaca mulatta',
            'sheep': 'Ovis aries',
            'southern pig-tailed macaque': 'Macaca nemestrina',
            'stump-tailed macaque': 'Macaca arctoides',
            'white-faced saki': 'Pithecia pithecia',
            'white-fronted spider monkey': 'Ateles belzebuth',
            'yellow baboon': 'Papio cynocephalus',
        }

    def reference_match_headers(self):
        return [
            'Reference Match',
        ]

    def get_mt_peptides(self):
        records = list(SeqIO.parse(self.input_fasta, "fasta"))
        if self.file_type == 'vcf':
            records_dict = {x.id.replace('MT.', ''): str(x.seq) for x in filter(lambda x: x.id.startswith('MT.'), records)}
        else:
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

    def extract_n_mer_from_fs(self, full_peptide, wt_peptide, epitope, subpeptide_position):
        #For frameshifts we want to test all downstream epitopes in the flanking region since they are all potentially novel
        flanking_sequence_length = self.match_length - 1
        start = subpeptide_position - 1 - flanking_sequence_length
        if start < 0:
            start = 0
        #This catches cases where the start position would cause too many leading wildtype amino acids, which would result
        #in false-positive reference matches
        if len(full_peptide) > len(wt_peptide):
            diff_position = [i for i in range(len(wt_peptide)) if wt_peptide[i] != full_peptide[i]][0]
        else:
            diff_position = [i for i in range(len(full_peptide)) if wt_peptide[i] != full_peptide[i]][0]
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
            processed_peptides = []
            reference_match_dict = defaultdict(list)
            for line in reader:
                if self.file_type == 'pVACbind':
                    epitope = line['Epitope Seq']
                    peptide = mt_records_dict[line['Mutation']]
                else:
                    epitope = line['MT Epitope Seq']
                    if self.file_type == 'vcf':
                        if line['Variant Type'] == 'FS':
                            peptide = self.extract_n_mer_from_fs(mt_records_dict[line['Index']], wt_records_dict[line['Index']], epitope, int(line['Sub-peptide Position']))
                        else:
                            mt_amino_acids = line['Mutation'].split('/')[1]
                            if mt_amino_acids == '-':
                                mt_amino_acids = ''
                            peptide = self.extract_n_mer(mt_records_dict[line['Index']], int(line['Sub-peptide Position']), int(line['Mutation Position']), len(mt_amino_acids))
                    else:
                        peptide = mt_records_dict[line['Index']]
                if peptide not in processed_peptides:
                    processed_peptides.append(peptide)
                    result_handle = NCBIWWW.qblast("blastp", "refseq_protein", peptide, entrez_query="{} [Organism]".format(self.species_to_organism[self.species]), word_size=min(self.match_length, 7), gapcosts='32767 32767')
                    for blast_record in NCBIXML.parse(result_handle):
                        if len(blast_record.alignments) > 0:
                            for alignment in blast_record.alignments:
                                for hsp in alignment.hsps:
                                    matches = re.split('\+| ', hsp.match)
                                    for match in matches:
                                        if len(match) >= self.match_length:
                                            reference_match_dict[peptide].append({
                                                'Hit ID': alignment.hit_id,
                                                'Hit Definition': alignment.hit_def,
                                                'Query Sequence': hsp.query,
                                                'Match Sequence': hsp.match,
                                                'Match Start': hsp.sbjct_start,
                                                'Match Stop': hsp.sbjct_end,
                                            })
                    sleep(10)
                if peptide in reference_match_dict:
                    line['Reference Match'] = True
                    metric_line = line.copy()
                    metric_line['Peptide'] = peptide
                    for alignment in reference_match_dict[peptide]:
                        metric_line.update(alignment)
                        metric_writer.writerow(metric_line)
                else:
                    line['Reference Match'] = False
                writer.writerow(line)
