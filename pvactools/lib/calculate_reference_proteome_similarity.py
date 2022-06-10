import csv
import argparse
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO, SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import shutil
import re
import os
import sys
from collections import defaultdict
from subprocess import run, DEVNULL, STDOUT
import tempfile
from time import sleep
import pymp

class CalculateReferenceProteomeSimilarity:
    '''
    Peforms blast search on the neoantigens found in the pipeline execution. 

    ...
    Parameters
    ----------
    input_file : str
        The path to <sample_name>.all_epitopes.tsv. It is the output of the MHC class I and II predictions. 

    input_fasta : str
        The path to a fasta file that contains the peptide sequences.

    output_file : str
        A path to a tsv file where the results of the BLAST calls are wrote to

    match_length : int
        The desired k-mer length

    species : str
        Species of organism blast is being performed on

    file_type : str
        The module you used to obtain MHC predictions:
            "pVACseq", "pVACbind" or "pVACfuse" 

    blastp_path : str
        A path to local install of the blastp tool

    blastp_db : str
        The name of the database to perform blast with

    n_threads : int
        The number of threads for multiprocessing


    Methods
    -------
    reference_match_headers()
        Returns the header for match result

    get_mt_peptides()
        Returns a dictionary of the mutant type peptides from the input fasta

    get_wt_peptides()
        Returns a dictionary of the wild type peptides from the input fasta

    extract_n_mer(full_peptide, subpeptide_position, mutation_position, mt_length)
        Returns n_mer from the full_peptide for not frameshift mutations

    extract_n_mer_from_fs(full_peptide, wt_peptide, epitope, subpeptide_position)
        returns n_mer from the full_peptide for frameshift mutations

    metric_headers()
        Returns the headers for the metric file

    _get_peptide(line, mt_records_dict, wt_records_dict)
        Returns the full_peptide and it's respective n_mer from the line in self.input_file

    _call_blast(self, full_peptide, peptide, p)
        Performs blast operations and returns tmp file that the results are written to.
        Note:
            TMP FILE NEEDS TO BE CLOSED OUTSIDE OF FUNCTION

    _needs_processing(full_peptide, processed_peptides)
        Returns true if protein has not been processed and false otherwise

    _generate_reference_match_dict(self, full_peptide, processed_peptides, peptide, p)
        Returns a dictionary that contains information about matches

    _write_outputs(input_fh, processed_peptides, mt_records_dict, wt_records_dict)
        Uses the blast records in processed_peptides to add results to information in the input_file and
        writes the new data to files.

    _get_unique_peptides(self, mt_records_dict, wt_records_dict)
        Creates a list of unique peptides from the input file 

    execute()
        Peforms the calculation of reference proteome similarity. The only method that should be 
        called from outside of the class
    '''
    def __init__(self, input_file, input_fasta, output_file, match_length=8, species='human', file_type='pVACseq', blastp_path=None, blastp_db='refseq_select_prot', n_threads=1):
        self.input_file = input_file
        self.input_fasta = input_fasta
        output_dir = os.path.dirname(output_file)
        os.makedirs(os.path.abspath(output_dir), exist_ok=True)
        self.output_file = output_file
        self.metric_file = "{}.reference_matches".format(output_file)
        self.match_length = match_length
        self.species = species
        self.n_threads = 1 if sys.platform == "darwin" else n_threads #pymp and requests not compatible on macOS 10.13+ for n_threads > 1
        self.file_type = file_type
        self.blastp_path = blastp_path
        self.blastp_db = blastp_db
        if self.blastp_db == 'refseq_select_prot' and self.species != 'human' and self.species != 'mouse':
            raise Exception("refseq_select_prot blastp database is only compatible with human and mouse species.")
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
        # Make a list of SeqRecords from the input_fasta
        records = list(SeqIO.parse(self.input_fasta, "fasta"))

        # Create mt record dictionary
        if self.file_type == 'pVACseq':
            records_dict = {x.id.replace('MT.', ''): str(x.seq) for x in filter(lambda x: x.id.startswith('MT.'), records)}
        else:
            records_dict = {x.id: str(x.seq) for x in records}
        return records_dict


    def get_wt_peptides(self):
        if self.file_type == 'pVACseq':
            # make a list of SeqRecords from the input_fasta
            records = list(SeqIO.parse(self.input_fasta, "fasta"))

            # Create wt record dictionary
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
        epitope_seq = 'MT Epitope Seq' if self.file_type == 'pVACseq' else 'Epitope Seq'
        return ['Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 'Transcript', epitope_seq, 'Peptide', 'Hit ID', 'Hit Definition', 'Query Sequence', 'Query Window', 'Match Sequence', 'Match Start', 'Match Stop']


    def _get_peptide(self, line, mt_records_dict, wt_records_dict):
        ## Get epitope, peptide and full_peptide
        if self.file_type == 'pVACbind' or self.file_type == 'pVACfuse':
            peptide = mt_records_dict[line['Mutation']]
            full_peptide = peptide
        else:
            epitope = line['MT Epitope Seq']
            full_peptide = mt_records_dict[line['Index']]

            # get peptide
            if line['Variant Type'] == 'FS':
                peptide = self.extract_n_mer_from_fs(full_peptide, wt_records_dict[line['Index']], epitope, int(line['Sub-peptide Position']))
            else:
                mt_amino_acids = line['Mutation'].split('/')[1]
                if mt_amino_acids == '-':
                    mt_amino_acids = ''
                peptide = self.extract_n_mer(full_peptide, int(line['Sub-peptide Position']), int(line['Mutation Position']), len(mt_amino_acids))
        return peptide, full_peptide


    def _call_blast(self, full_peptide, p):
        if self.blastp_path is not None: # if blastp installed locally, perform BLAST with it

            # create a SeqRecord of full_peptide and write it to a tmp file
            record = SeqRecord(Seq(full_peptide, IUPAC.protein), id="1", description="")
            tmp_peptide_fh = tempfile.NamedTemporaryFile('w', delete=False)
            SeqIO.write([record], tmp_peptide_fh.name, "fasta")

            # configure args for local blastp, run it and put results in new tmp file
            arguments = [self.blastp_path, '-query', tmp_peptide_fh.name, '-db', self.blastp_db, '-outfmt', '16', '-word_size', str(min(self.match_length, 7)), '-gapopen', '32767', '-gapextend', '32767']
            result_handle = tempfile.NamedTemporaryFile(delete=False)
            response = run(arguments, stdout=result_handle, check=True)
            result_handle.seek(0)
            tmp_peptide_fh.close()

        else: # else perform BLAST with api
            with p.lock: # stagger calls to qblast
                if not os.environ.get('TEST_FLAG') or os.environ.get('TEST_FLAG') == '0': # we don't need to sleep during testing since this is mocked and not actually calling the API
                    sleep(10)
            result_handle = NCBIWWW.qblast("blastp", self.blastp_db, full_peptide, entrez_query="{} [Organism]".format(self.species_to_organism[self.species]), word_size=min(self.match_length, 7), gapcosts='32767 32767')

        return result_handle


    def _generate_reference_match_dict(self, blast_records, peptide):

        reference_match_dict = defaultdict(list)
        for blast_record in blast_records:
            if len(blast_record.alignments) > 0: # if there is at least one alignment
                for alignment in blast_record.alignments:
                    if alignment.title.endswith(" [{}]".format(self.species_to_organism[self.species])):
                        for hsp in alignment.hsps: # High-scoring Segment Pair
                            # create list of strings that represent matching windows
                            matches = re.split('\+| ', hsp.match)
                            for match in matches:
                                # 'windows' of query peptides that match subject peptides
                                windows = [match[i:i+self.match_length] for i in range(len(match)-(self.match_length-1))]
                                for window in windows: 
                                    if window in peptide:
                                        reference_match_dict[peptide].append({
                                            'Hit ID': alignment.hit_id,
                                            'Hit Definition': alignment.hit_def,
                                            'Query Sequence': hsp.query,
                                            'Query Window'  : window,
                                            'Match Sequence': hsp.match,
                                            'Match Start': hsp.sbjct_start,
                                            'Match Stop': hsp.sbjct_end,
                                        })
        return reference_match_dict


    def _write_outputs(self, processed_peptides, mt_records_dict, wt_records_dict):

        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh, open(self.metric_file, 'w') as metric_fh:
            reader = csv.DictReader(input_fh, delimiter="\t")
            writer = csv.DictWriter(output_fh, delimiter="\t", fieldnames=reader.fieldnames + self.reference_match_headers(), extrasaction='ignore')
            metric_writer = csv.DictWriter(metric_fh, delimiter="\t", fieldnames=self.metric_headers(), extrasaction='ignore')
            writer.writeheader()
            metric_writer.writeheader()

            for line in reader:
                peptide, full_peptide = self._get_peptide(line, mt_records_dict, wt_records_dict)

                blast_records = processed_peptides[full_peptide]
                reference_match_dict = self._generate_reference_match_dict(blast_records, peptide)

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

    def _get_unique_peptides(self, mt_records_dict, wt_records_dict):
        unique_peptides = set()

        with open(self.input_file) as input_fh:
            reader = csv.DictReader(input_fh, delimiter='\t')
            for line in reader:
                _, full_peptide = self._get_peptide(line, mt_records_dict, wt_records_dict)
                unique_peptides.add(full_peptide)

        return list(unique_peptides)


    def execute(self):

        if self.species not in self.species_to_organism:
            print("Species {} not supported for Reference Proteome Similarity search. Skipping.".format(self.species))
            shutil.copy(self.input_file, self.output_file)
            return

        mt_records_dict = self.get_mt_peptides()
        wt_records_dict = self.get_wt_peptides()

        unique_peptides = pymp.shared.list(self._get_unique_peptides(mt_records_dict, wt_records_dict))
        processed_peptides = pymp.shared.dict()

        with pymp.Parallel(self.n_threads) as p:
            for i in p.range(len(unique_peptides)):

                full_peptide = unique_peptides[i]

                result_handle = self._call_blast(full_peptide, p)

                blast_records = [x for x in NCBIXML.parse(result_handle)]

                with p.lock:
                    processed_peptides[full_peptide] = blast_records
                result_handle.close()

        self._write_outputs(processed_peptides, mt_records_dict, wt_records_dict)

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s calculate_reference_proteome_similarity' % tool,
            description="Blast peptides against the reference proteome.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="Input filtered or all_epitopes file with predicted epitopes."
        )
        parser.add_argument(
            'input_fasta',
            help="For pVACbind, the original input FASTA file. "
            + "For pVACseq and pVACfuse a FASTA file with mutant peptide sequences for each variant isoform. "
            + "This file can be found in the same directory as the input filtered/all_epitopes file. "
            + "Can also be generated by running `pvacseq|pvacfuse generate_protein_fasta`.")
        parser.add_argument(
            'output_file',
            help="Output TSV filename for putative neoepitopes."
        )
        parser.add_argument(
            '--match-length',
            default=8,
            type=int,
            help="The desired matching epitope length."
        )
        parser.add_argument(
            '--species',
            default='human',
            help="The species of the input file."
        )
        parser.add_argument(
            '--blastp-path',
            default=None,
            help="Blastp installation path.",
        )
        parser.add_argument(
            '--blastp-db',
            choices=['refseq_select_prot', 'refseq_protein'],
            default='refseq_select_prot',
            help="The blastp database to use.",
        )
        parser.add_argument(
            "-t", "--n-threads",type=int,
            default=1,
            help="Number of threads to use for parallelizing BLAST calls.",
        )
        return parser
