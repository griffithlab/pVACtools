import csv
import argparse
import gzip
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
from itertools import groupby
import json

from pvactools.lib.run_utils import *

class CalculateReferenceProteomeSimilarity:
    '''
    Peforms blast search on the neoantigens found in the pipeline execution. 

    ...
    Parameters
    ----------
    input_file : str
        The path to <sample_name>.all_epitopes.tsv, <sample_name>.filtered.tsv, or <sample_name>.all_epitopes.aggregated.tsv. It is the output of the MHC class I and II predictions. 

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

    peptide_fasta: str
        The path to a FASTA file with transcript peptide sequences to search for matches

    n_threads : int
        The number of threads for multiprocessing


    Methods
    -------
    reference_match_headers(fieldnames)
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

    _match_from_peptide_fasta(self, full_peptide)
        Finds matches for the windows in full_peptide in the peptide fasta

    _needs_processing(full_peptide, processed_peptides)
        Returns true if protein has not been processed and false otherwise

    _generate_reference_match_dict_from_blast_records(self, blast_records, peptide)
        Returns a dictionary that contains information about matches obtained from blast

    _generate_reference_match_dict_from_peptide_fasta_results(self, results, peptide, transcript)
        Returns a dictionary that contains information about matches obtained from the peptide fasta

    _write_outputs(input_fh, processed_peptides, mt_records_dict, wt_records_dict)
        Uses the blast records in processed_peptides to add results to information in the input_file and
        writes the new data to files.

    _get_unique_peptides(self, mt_records_dict, wt_records_dict)
        Creates a list of unique peptides from the input file 

    execute()
        Peforms the calculation of reference proteome similarity. The only method that should be 
        called from outside of the class
    '''
    def __init__(self, input_file, input_fasta, output_file, match_length=8, species='human', file_type='pVACseq', blastp_path=None, blastp_db='refseq_select_prot', peptide_fasta=None, n_threads=1, aggregate_metrics_file=None):
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
        self.peptide_fasta = peptide_fasta
        self.aggregate_metrics_file = aggregate_metrics_file
        if self.aggregate_metrics_file:
            (dir_name, full_file_name) = os.path.split(output_file)
            (file_name, extension) = os.path.splitext(full_file_name)
            if extension != '.tsv':
                raise Exception("Output file name is expected to be a .tsv file.")
            self.output_aggregate_metrics_file = output_file.replace('.tsv', '.metrics.json')
            with open(self.aggregate_metrics_file, 'r') as fh:
                self.aggregate_metrics = json.loads(fh.read())
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


    def reference_match_headers(self, fieldnames):
        if self._input_tsv_type(fieldnames) == 'aggregated':
            if 'Ref Match' not in fieldnames:
                fieldnames.insert(len(fieldnames)-1, 'Ref Match')
        else:
            if 'Reference Match' not in fieldnames:
                fieldnames.insert(len(fieldnames), 'Reference Match')
        return fieldnames


    def get_mt_peptides(self):
        # Make a list of SeqRecords from the input_fasta
        records = list(SeqIO.parse(self.input_fasta, "fasta"))

        # Create mt record dictionary
        if self.file_type == 'pVACseq':
            records_dict = {re.sub('^%s' % "MT\.", "", x.id): str(x.seq) for x in filter(lambda x: x.id.startswith('MT.'), records)}
        else:
            records_dict = {x.id: str(x.seq) for x in records}
        return records_dict


    def get_wt_peptides(self):
        if self.file_type == 'pVACseq':
            # make a list of SeqRecords from the input_fasta
            records = list(SeqIO.parse(self.input_fasta, "fasta"))

            # Create wt record dictionary
            records_dict = {re.sub('^%s' % "WT\.", "", x.id): str(x.seq) for x in filter(lambda x: x.id.startswith('WT.'), records)}
        else:
            return {}
        return records_dict


    def extract_n_mer(self, full_peptide, subpeptide_position, mutation_position, mt_length, variant_type):
        #For non-frameshifts this ensures that we only test match_length epitopes that overlap the mutation
        #If we extract a larger region, we will get false-positive matches against the reference proteome
        #from the native wildtype portion of the peptide
        flanking_sequence_length = self.match_length - 1
        if variant_type == 'inframe_del':
            mt_start = subpeptide_position + (mutation_position)
        else:
            mt_start = subpeptide_position + (mutation_position-1)
        start = mt_start - flanking_sequence_length
        if start < 0:
            start = 0
        if variant_type == 'inframe_del':
            end = mt_start + flanking_sequence_length
        else:
            end = mt_start + mt_length + flanking_sequence_length
        return full_peptide[start:end]


    def extract_n_mer_from_fs(self, full_peptide, wt_peptide, epitope, subpeptide_position):
        #For frameshifts we want to test all downstream epitopes in the flanking region since they are all potentially novel
        flanking_sequence_length = self.match_length - 1
        start = subpeptide_position - flanking_sequence_length
        if start < 0:
            start = 0
        #This catches cases where the start position would cause too many leading wildtype amino acids, which would result
        #in false-positive reference matches
        if len(full_peptide) > len(wt_peptide):
            diffs = [i for i in range(len(wt_peptide)) if wt_peptide[i] != full_peptide[i]]
            if diffs == []:
                diffs = [len(wt_peptide)]
        else:
            diffs = [i for i in range(len(full_peptide)) if wt_peptide[i] != full_peptide[i]]
            if diffs == []:
                diffs = [len(full_peptide)]
        diff_position = diffs[0]
        min_start = diff_position - self.match_length + 1 
        if min_start > start:
            start = min_start
        end = start + flanking_sequence_length + len(epitope) + flanking_sequence_length
        return full_peptide[start:end]


    def metric_headers(self):
        epitope_seq = 'MT Epitope Seq' if self.file_type == 'pVACseq' else 'Epitope Seq'
        if self.file_type == 'pVACseq':
            return ['Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 'Transcript', epitope_seq, 'Peptide', 'Hit ID', 'Hit Definition', 'Match Window', 'Match Sequence', 'Match Start', 'Match Stop']
        else:
            return ['ID', epitope_seq, 'Peptide', 'Hit ID', 'Hit Definition', 'Match Window', 'Match Sequence', 'Match Start', 'Match Stop']


    def _input_tsv_type(self, line):
        if self.file_type == 'pVACseq' and 'Index' in line:
            return 'full'
        elif self.file_type != 'pVACseq' and 'Mutation' in line:
            return 'full'
        else:
            return 'aggregated'

    def _get_full_peptide(self, line, mt_records_dict, wt_records_dict):
        for record_id in mt_records_dict.keys():
            (rest_record_id, variant_type, aa_change) = record_id.rsplit(".", 2)
            (count, gene, transcript) = rest_record_id.split(".", 2)
            (parsed_aa_change, pos, wt_aa, mt_aa) = index_to_aggregate_report_aa_change(aa_change, variant_type)
            if line['Best Transcript'] == transcript and line['AA Change'] == parsed_aa_change:
                return (mt_records_dict[record_id], wt_records_dict[record_id], variant_type, mt_aa, wt_aa)
        print("Unable to find full_peptide for variant {}".format(line['ID']))
        return (None, None, variant_type, mt_aa, wt_aa)

    def _get_peptide(self, line, mt_records_dict, wt_records_dict):
        ## Get epitope, peptide and full_peptide
        if self.file_type == 'pVACbind' or self.file_type == 'pVACfuse':
            if self._input_tsv_type(line) == 'aggregated':
                peptide = mt_records_dict[line['ID']]
            else:
                peptide = mt_records_dict[line['Mutation']]
            full_peptide = peptide
        else:
            if self._input_tsv_type(line) == 'aggregated':
                epitope = line['Best Peptide']
                (full_peptide, wt_peptide, variant_type, mt_amino_acids, wt_amino_acids) = self._get_full_peptide(line, mt_records_dict, wt_records_dict)
                if full_peptide is None:
                    return None, None
                if variant_type != 'FS':
                    if line['Pos'] == 'NA':
                        mt_pos = None
                        for i,(wt_aa,mt_aa) in enumerate(zip(wt_peptide,full_peptide)):
                            if wt_aa != mt_aa:
                                mt_pos = i
                                break
                        if mt_pos is None:
                            return None, full_peptide
                    else:
                        mt_pos = int(line['Pos'].split('-')[0])
            else:
                epitope = line['MT Epitope Seq']
                full_peptide = mt_records_dict[line['Index']]
                wt_peptide = wt_records_dict[line['Index']]
                variant_type = line['Variant Type']
                if variant_type != 'FS':
                    mt_pos = int(line['Mutation Position'].split('-')[0])
                    (wt_amino_acids, mt_amino_acids) = line['Mutation'].split('/')

            # get peptide
            subpeptide_position = full_peptide.index(epitope)
            if variant_type == 'FS':
                peptide = self.extract_n_mer_from_fs(full_peptide, wt_peptide, epitope, subpeptide_position)
            else:
                if mt_amino_acids == '-':
                    mt_amino_acids = ''
                if len(mt_amino_acids) == len(wt_amino_acids) and len(mt_amino_acids) > 1:
                    #remove leading and trailing identical amino acids
                    shortened_mt_amino_acids = ""
                    for mt_aa, wt_aa in zip(mt_amino_acids, wt_amino_acids):
                        if wt_aa != mt_aa:
                            shortened_mt_amino_acids += mt_aa
                else:
                    shortened_mt_amino_acids = mt_amino_acids
                peptide = self.extract_n_mer(full_peptide, subpeptide_position, mt_pos, len(shortened_mt_amino_acids), variant_type)
        return peptide, full_peptide


    def _call_blast(self, full_peptide, p):
        word_size = min(self.match_length, 7, int(len(full_peptide)/2))
        if self.blastp_path is not None: # if blastp installed locally, perform BLAST with it

            # create a SeqRecord of full_peptide and write it to a tmp file
            record = SeqRecord(Seq(full_peptide, IUPAC.protein), id="1", description="")
            tmp_peptide_fh = tempfile.NamedTemporaryFile('w', delete=False)
            SeqIO.write([record], tmp_peptide_fh.name, "fasta")

            # configure args for local blastp, run it and put results in new tmp file
            arguments = [self.blastp_path, '-query', tmp_peptide_fh.name, '-db', self.blastp_db, '-outfmt', '16', '-word_size', str(word_size), '-gapopen', '32767', '-gapextend', '32767']
            result_handle = tempfile.NamedTemporaryFile(delete=False)
            response = run(arguments, stdout=result_handle, check=True)
            result_handle.seek(0)
            tmp_peptide_fh.close()

        else: # else perform BLAST with api
            with p.lock: # stagger calls to qblast
                if not os.environ.get('TEST_FLAG') or os.environ.get('TEST_FLAG') == '0': # we don't need to sleep during testing since this is mocked and not actually calling the API
                    sleep(10)
            result_handle = NCBIWWW.qblast("blastp", self.blastp_db, full_peptide, entrez_query="{} [Organism]".format(self.species_to_organism[self.species]), word_size=word_size, gapcosts='32767 32767', hitlist_size=500)

        return result_handle


    def _match_from_peptide_fasta(self, full_peptide):
        results = []
        with gzip.open(self.peptide_fasta, 'rt') as handle:
            transcript_sequences = SeqIO.parse(handle, "fasta")
            for transcript_seq in transcript_sequences:
                seq = str(transcript_seq.seq)
                for i in range(0, len(full_peptide)-self.match_length+1):
                    epitope = full_peptide[i:i+self.match_length]
                    if epitope in seq:
                        results.append([transcript_seq, epitope, seq.index(epitope)])
        return results


    def _generate_reference_match_dict_from_blast_records(self, blast_records, peptide):
        reference_match_dict = []
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
                                        reference_match_dict.append({
                                            'Hit ID': alignment.hit_id,
                                            'Hit Definition': alignment.hit_def,
                                            'Match Window'  : window,
                                            'Match Sequence': match,
                                            'Match Start': match.index(window) + 1 ,
                                            'Match Stop': match.index(window) + 1 + len(window),
                                        })
        return self._combine_reference_match_entries(reference_match_dict)


    def _generate_reference_match_dict_from_peptide_fasta_results_for_pvacbind(self, results, peptide):
        reference_match_dict = []
        for transcript_seq, epitope, match_start in results:
            reference_match_dict.append({
                'Hit ID': transcript_seq.id,
                'Hit Definition': transcript_seq.description,
                'Match Window'  : epitope,
                'Match Sequence': str(transcript_seq.seq),
                'Match Start': match_start + 1,
                'Match Stop': match_start + len(epitope) + 1,
            })
        return self._combine_reference_match_entries(reference_match_dict)

    def _generate_reference_match_dict_from_peptide_fasta_results(self, results, peptide, transcript):
        reference_match_dict = []
        for transcript_seq, epitope, match_start in results:
            if transcript in transcript_seq.description:
                continue
            reference_match_dict.append({
                'Transcript': transcript,
                'Hit ID': transcript_seq.id,
                'Hit Definition': transcript_seq.description,
                'Match Window'  : epitope,
                'Match Sequence': str(transcript_seq.seq),
                'Match Start': match_start + 1,
                'Match Stop': match_start + len(epitope) + 1,
            })
        return self._combine_reference_match_entries(reference_match_dict)

    def _combine_reference_match_entries(self, reference_match_dict):
        combined_matches = []
        for hit_id, hit_reference_matches in groupby(reference_match_dict,key=lambda x:x['Hit ID']):
            hit_reference_matches = sorted(list(hit_reference_matches), key=lambda d: d['Match Start'])
            first_match = None
            for index, match in enumerate(hit_reference_matches):
                if first_match is None:
                    combined_match = match.copy()
                    match_start = combined_match['Match Start']
                    first_match = match_start
                if index == len(hit_reference_matches) - 1:
                    match_stop = hit_reference_matches[index]['Match Stop']
                    combined_match['Match Stop'] = match_stop
                    combined_match['Match Window'] = match['Match Sequence'][match_start-1:match_stop-1]
                    combined_matches.append(combined_match)
                elif hit_reference_matches[index+1]['Match Start'] != match['Match Start'] + 1:
                    match_stop = match['Match Stop']
                    combined_match['Match Stop'] = match_stop
                    combined_match['Match Window'] = match['Match Sequence'][match_start-1:match_stop-1]
                    combined_matches.append(combined_match)
                    first_match = None
        unique_combined_matches = [dict(s) for s in set(frozenset(d.items()) for d in combined_matches)]
        combined_reference_matches = sorted(unique_combined_matches, key=lambda x: (x['Match Stop'] - x['Match Start'], x['Hit ID'], x['Match Start']), reverse=True)
        return combined_reference_matches


    def _write_outputs(self, processed_peptides, mt_records_dict, wt_records_dict):

        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh, open(self.metric_file, 'w') as metric_fh:
            reader = csv.DictReader(input_fh, delimiter="\t")
            writer = csv.DictWriter(output_fh, delimiter="\t", fieldnames=self.reference_match_headers(reader.fieldnames.copy()), extrasaction='ignore')
            metric_writer = csv.DictWriter(metric_fh, delimiter="\t", fieldnames=self.metric_headers(), extrasaction='ignore')
            writer.writeheader()
            metric_writer.writeheader()

            for line in reader:
                peptide, full_peptide = self._get_peptide(line, mt_records_dict, wt_records_dict)

                if self.peptide_fasta:
                    if peptide is None:
                        if self._input_tsv_type(line) == 'aggregated':
                            line['Ref Match'] = 'Not Run'
                            if self.aggregate_metrics_file:
                                self.aggregate_metrics[line['ID']]['reference_matches'] = {
                                    'count': 0,
                                    'query_peptide': peptide,
                                    'matches': []
                                }
                        else:
                            line['Reference Match'] = 'Not Run'
                        writer.writerow(line)
                        continue

                    results = processed_peptides[peptide]
                else:
                    results = processed_peptides[full_peptide]

                if self.peptide_fasta:
                    if self.file_type == 'pVACseq' or self.file_type == 'pVACfuse':
                        if self._input_tsv_type(line) == 'full':
                            transcript = line['Transcript']
                        else:
                            transcript = line['Best Transcript']
                        reference_matches = self._generate_reference_match_dict_from_peptide_fasta_results(results, peptide, transcript)
                    else:
                        reference_matches = self._generate_reference_match_dict_from_peptide_fasta_results_for_pvacbind(results, peptide)
                else:
                    reference_matches = self._generate_reference_match_dict_from_blast_records(results, peptide)

                if len(reference_matches) > 0:
                    if self._input_tsv_type(line) == 'aggregated':
                        line['Ref Match'] = True
                    else:
                        line['Reference Match'] = True
                    metric_line = line.copy()
                    if self._input_tsv_type(line) == 'aggregated' and self.file_type == 'pVACseq':
                        (chromosome, start, stop, ref, alt) = line['ID'].split('-')
                        metric_line['Chromosome'] = chromosome
                        metric_line['Start'] = start
                        metric_line['Stop'] = stop
                        metric_line['Reference'] = ref
                        metric_line['Variant'] = alt
                    metric_line['Peptide'] = peptide
                    if self._input_tsv_type(line) == 'aggregated':
                        epitope_seq_header = 'MT Epitope Seq' if self.file_type == 'pVACseq' else 'Epitope Seq'
                        metric_line[epitope_seq_header] = line['Best Peptide']
                    metric_lines = []
                    for alignment in reference_matches:
                        metric_line_alignment = metric_line.copy()
                        metric_line_alignment.update(alignment)
                        metric_lines.append(metric_line_alignment)
                    if self.aggregate_metrics_file:
                        matches = []
                        for query_window, hit_reference_matches in groupby(metric_lines,key=lambda x:x['Match Window']):
                            hit_reference_matches = list(hit_reference_matches)
                            gene_regex = '^.*gene_symbol:([0-9|A-Z]+).*$'
                            transcript_regex = '^.*transcript:(ENS[0-9|A-Z|.]+).*$'
                            gene_p = re.compile(gene_regex)
                            transcript_p = re.compile(transcript_regex)
                            genes = []
                            transcripts = []
                            for definition in [l['Hit Definition'] for l in hit_reference_matches]:
                                m = gene_p.match(definition)
                                if m:
                                    genes.append(m.group(1))
                                m = transcript_p.match(definition)
                                if m:
                                    transcripts.append(m.group(1))
                            if len(genes) > 0:
                                matches.append({
                                    'Matched Peptide': query_window,
                                    'Genes': ", ".join(sorted(list(set(genes)))),
                                    'Transcripts': ", ".join(sorted(list(set(transcripts)))),
                                    'Hit IDs': [l['Hit ID'] for l in hit_reference_matches],
                                })
                            else:
                                matches.append({
                                    'Matched Peptide': query_window,
                                    'Hit IDs': [l['Hit ID'] for l in hit_reference_matches],
                                    'Hit Definitions': [l['Hit Definition'] for l in hit_reference_matches],
                                })
                        m = {
                            'count': len(metric_lines),
                            'query_peptide': metric_lines[0]['Peptide'],
                            'matches': matches
                        }
                        self.aggregate_metrics[metric_lines[0]['ID']]['reference_matches'] = m
                    metric_writer.writerows(metric_lines)

                else:
                    if self._input_tsv_type(line) == 'aggregated':
                        line['Ref Match'] = False
                        if self.aggregate_metrics_file:
                            self.aggregate_metrics[line['ID']]['reference_matches'] = {
                                'count': 0,
                                'query_peptide': peptide,
                                'matches': []
                            }
                    else:
                        line['Reference Match'] = False
                writer.writerow(line)
        if self.aggregate_metrics_file:
            with open(self.output_aggregate_metrics_file, 'w') as fh:
                json.dump(self.aggregate_metrics, fh, indent=2, separators=(',', ': '))

    def _get_unique_peptides(self, mt_records_dict, wt_records_dict):
        unique_peptides = set()

        with open(self.input_file) as input_fh:
            reader = csv.DictReader(input_fh, delimiter='\t')
            for line in reader:
                peptide, full_peptide = self._get_peptide(line, mt_records_dict, wt_records_dict)
                if self.peptide_fasta:
                    if peptide is not None:
                        unique_peptides.add(peptide)
                else:
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

                if self.peptide_fasta:
                    results = self._match_from_peptide_fasta(full_peptide)
                else:
                    result_handle = self._call_blast(full_peptide, p)
                    results = [x for x in NCBIXML.parse(result_handle)]
                    result_handle.close()

                with p.lock:
                    processed_peptides[full_peptide] = results

        self._write_outputs(processed_peptides, mt_records_dict, wt_records_dict)

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s calculate_reference_proteome_similarity' % tool,
            description="Identify which epitopes in a pVACseq|pVACfuse|pVACbind report file have matches in the reference proteome using either BLASTp or a checking directly against a reference proteome FASTA.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="Input filtered, all_epitopes, or aggregated report file with predicted epitopes."
        )
        parser.add_argument(
            'input_fasta',
            help="For pVACbind, the original input FASTA file. "
            + "For pVACseq and pVACfuse a FASTA file with mutant peptide sequences for each variant isoform. "
            + "This file can be found in the same directory as the input filtered.tsv/all_epitopes.tsv file. "
            + "Can also be generated by running `pvacseq|pvacfuse generate_protein_fasta`.")
        parser.add_argument(
            'output_file',
            help="Output TSV filename of report file with epitopes with reference matches marked."
        )
        parser.add_argument(
            '--match-length',
            default=8,
            type=int,
            help="The minimum number of consecutive amino acids that need to match."
        )
        parser.add_argument(
            '--species',
            default='human',
            help="The species of the data in the input file."
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
            '--peptide-fasta',
            help="A reference peptide FASTA file to use for finding reference matches instead of blastp."
        )
        parser.add_argument(
            "-t", "--n-threads",type=int,
            default=1,
            help="Number of threads to use for parallelizing BLAST calls.",
        )
        if tool == 'pvacseq':
            parser.add_argument(
                "-m", "--aggregate-metrics-file",
                help="When running with the aggregate report as an input tsv, optionally provide the metrics.json "
                    + "file to update with detailed reference match data for display in pVACview."
            )
        return parser
