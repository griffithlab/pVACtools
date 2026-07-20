import tempfile
import os
import shutil
import csv
import re
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pvactools.lib.calculate_manufacturability import CalculateManufacturability
from pvactools.lib.run_utils import *

class GenerateProteinFasta:
    def __init__(self, **kwargs):
        self.sample_name = kwargs.pop('sample_name', 'tmp')
        if self.sample_name is None:
            self.sample_name = 'tmp'
        self.pass_only = kwargs.pop('pass_only', False)
        self.biotypes = kwargs.pop('biotypes', ['protein_coding'])
        self.allow_incomplete_transcripts = kwargs.pop('allow_incomplete_transcripts', False)
        self.flanking_sequence_length = kwargs.pop('flanking_sequence_length')
        self.mutant_only = kwargs.pop('mutant_only', False)
        self.aggregate_report_evaluation = kwargs.pop('aggregate_report_evaluation', ['Accept'])
        self.temp_dir = tempfile.mkdtemp()
        self.fasta_file_path = kwargs.pop('fasta_file_path', os.path.join(self.temp_dir, f"{self.sample_name}.transcripts.fa"))
        self.trimmed_fasta_file_path = kwargs.pop('trimmed_fasta_file_path', os.path.join(self.temp_dir, f"{self.sample_name}.transcripts.trimmed.fa"))
        self.filtered_fasta_file_path = os.path.join(self.temp_dir, f"{self.sample_name}.transcripts.filtered.fa")
        self.input_tsv = kwargs.pop('input_tsv', None)
        self.output_file = kwargs.pop('output_file', None)

    def execute(self):
        self.generate_fasta()
        self.trim_sequences()
        self.filter_fasta()
        shutil.copy(self.filtered_fasta_file_path, self.output_file)
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        manufacturability_file = "{}.manufacturability.tsv".format(self.output_file)
        print("Calculating Manufacturability Metrics")
        CalculateManufacturability(self.output_file, manufacturability_file, 'fasta').execute()
        print("Completed")

    def generate_fasta(self):
        raise Exception("Implement in child class")

    def trim_sequences(self):
        raise Exception("Implement in child class")
    def parse_input_tsv(self):
        if self.input_tsv is None:
            return (None, None)
        indexes = []
        with open(self.input_tsv, 'r') as fh:
            reader = csv.DictReader(fh, delimiter = "\t")
            if 'Best Peptide' in reader.fieldnames:
                for line in reader:
                    indexes.append(line)
                file_type = 'aggregated'
            else:
                for line in reader:
                    indexes.append(line['Index'])
                file_type = 'full'
        return (indexes, file_type)

    def filter_fasta(self):
        if self.input_tsv is None:
            shutil.copy(self.trimmed_fasta_file_path, self.filtered_fasta_file_path)
        else:
            print("Filtering Variant Peptide FASTA")
            (tsv_indexes, file_type) = self.parse_input_tsv()

            output_records = []
            for record in SeqIO.parse(self.trimmed_fasta_file_path, "fasta"):
                record_id = record.id.split('.', 1)[1]
                description = ""
                if file_type == 'full':
                    if record_id not in tsv_indexes:
                        continue
                else:
                    matches = [r for r in tsv_indexes if r['Index'] == record_id and r['Evaluation'] in self.aggregate_report_evaluation]
                    if len(matches) == 0:
                        continue
                    elif len(matches) > 0 and (record.id.startswith('MT.') or record.id.startswith('ALT.')):
                        description = json.dumps({ 'Best Peptide': matches[0]['Best Peptide'] })
                new_record = SeqRecord(record.seq, id=record.id, description=description)
                output_records.append(new_record)

            ordered_output_records = []
            for tsv_index in tsv_indexes:
                if file_type == 'full':
                    records = [r for r in output_records if r.id.split('.', 1)[1] == tsv_index]
                else:
                    records = [r for r in output_records if r.id.split('.', 1)[1] == tsv_index['Index']]
                ordered_output_records.extend(records)
            output_records = ordered_output_records

            SeqIO.write(output_records, self.filtered_fasta_file_path, "fasta")
            print("Completed")

class PvacseqGenerateProteinFasta(GenerateProteinFasta):
    def __init__(self, **kwargs):
        self.input_vcf = kwargs.pop('input_vcf', None)
        self.phased_proximal_variants_vcf = kwargs.pop('phased_proximal_variants_vcf', None)
        self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', 1000)
        super().__init__(**kwargs)

    def generate_fasta(self):
        from pvactools.lib.variant_to_kmer_pipeline import VariantToKmerPipeline
        params = {
            'output_dir'                  : self.temp_dir,
            'input_file'                  : self.input_vcf,
            'sample_name'                 : self.sample_name,
            'pass_only'                   : self.pass_only,
            'proximal_variants_vcf'       : self.phased_proximal_variants_vcf,
            'biotypes'                    : self.biotypes,
            'allow_incomplete_transcripts': self.allow_incomplete_transcripts,
            'downstream_sequence_length'  : self.downstream_sequence_length,
            'flanking_bases'              : self.flanking_sequence_length,
        }
        pipeline = VariantToKmerPipeline(**params)
        pipeline.generate_fasta()

    def trim_sequences(self):
        print("Trimming Variant Peptide FASTA")
        records = {}
        keys = set()
        for record in SeqIO.parse(self.fasta_file_path, "fasta"):
            records[record.id] = str(record.seq)
            keys.add(record.id.split('.', 1)[1])

        output_records = []
        for key in sorted(keys, key=lambda x: int(x.split('.', 1)[0])):
            mt_seq = records[f"MT.{key}"]
            wt_seq = records[f"WT.{key}"]
            _, variant_type, aa_change = key.rsplit('.', 2)
            position = int(re.split('[A-Z|-]', aa_change)[0])
            start_position = position - self.flanking_sequence_length - 1
            end_position = position + self.flanking_sequence_length
            if variant_type == 'missense':
                if start_position < 0:
                    start_position = 0
                trimmed_mt_seq = mt_seq[start_position:end_position]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    trimmed_wt_seq = wt_seq[start_position:end_position]
                    output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))
            elif variant_type == 'FS':
                if start_position < 0:
                    start_position = 0
                trimmed_mt_seq = mt_seq[start_position:]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    trimmed_wt_seq = wt_seq[start_position:end_position]
                    output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))
            elif variant_type == 'inframe_del':
                match = re.match(r"\d+(?:-\d+)?([A-Z]+)/([A-Z]+|-)", aa_change)
                wt_aa, mt_aa = match.groups()
                if wt_aa.startswith(mt_aa):
                    start_position = start_position + 1
                else:
                    end_position = end_position - 1
                if start_position < 0:
                    start_position = 0
                trimmed_mt_seq = mt_seq[start_position:(end_position)]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    offset = len(wt_seq) - len(mt_seq)
                    trimmed_wt_seq = wt_seq[start_position:(end_position + offset)]
                    output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))
            else:
                match = re.match(r"\d+(?:-\d+)?([A-Z]+|-)/([A-Z]+)", aa_change)
                wt_aa, mt_aa = match.groups()
                if mt_aa.startswith(wt_aa):
                    start_position = start_position + 1
                else:
                    end_position = end_position - 1
                if start_position < 0:
                    start_position = 0
                offset = len(mt_seq) - len(wt_seq)
                trimmed_mt_seq = mt_seq[start_position:(end_position + offset)]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    trimmed_wt_seq = wt_seq[start_position:(end_position)]
                    output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))

        SeqIO.write(output_records, self.trimmed_fasta_file_path, "fasta")
        print("Completed")

class PvacspliceGenerateProteinFasta(GenerateProteinFasta):
    def __init__(self, **kwargs):
        self.input_file = kwargs.pop('input_file', None)
        self.annotated_vcf = kwargs.pop('annotated_vcf', None)
        self.ref_fasta = kwargs.pop('ref_fasta', None)
        self.gtf_file = kwargs.pop('gtf_file', None)
        self.junction_score = kwargs.pop('junction_score', 10)
        self.variant_distance = kwargs.pop('variant_distance', 100)
        self.anchor_types = kwargs.pop('anchor_types', ['A', 'D', 'NDA'])
        #self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', 1000)
        super().__init__(**kwargs)

    def generate_fasta(self):
        from pvactools.lib.junction_to_kmer_pipeline import JunctionToKmerPipeline
        junction_arguments = {
            'input_file_type'                  : 'junctions',
            'junctions_dir'                    : self.temp_dir,
            'input_file'                       : self.input_file,
            'gtf_file'                         : self.gtf_file,
            'save_gtf'                         : False,
            'sample_name'                      : self.sample_name,
            'ref_fasta'                        : self.ref_fasta,
            'annotated_vcf'                    : self.annotated_vcf,
            'pass_only'                        : self.pass_only,
            'biotypes'                         : self.biotypes,
            'allow_incomplete_transcripts'     : self.allow_incomplete_transcripts,
            'junction_score'                   : self.junction_score,
            'variant_distance'                 : self.variant_distance,
            'anchor_types'                     : self.anchor_types,
            'normal_sample_name'               : None,
            'keep_tmp_files'                   : False,
            'class_i_epitope_length'           : [],
            'class_ii_epitope_length'          : [],
            'class_i_hla'                      : [],
            'class_ii_hla'                     : [],
        }

        pipeline = JunctionToKmerPipeline(**junction_arguments)
        pipeline.generate_fasta()

        return pipeline.create_file_path('fasta')

    def trim_sequences(self):
        records = {}
        keys = set()
        for record in SeqIO.parse(self.fasta_file_path, "fasta"):
            records[record.id] = record.seq
            keys.add(record.id.split('.', 1)[1])

        output_records = []
        for key in sorted(keys, key=lambda x: int(x.split('.', 1)[0])):
            mt_sequence = records[f"ALT.{key}"]
            wt_sequence = records[f"WT.{key}"]
            if mt_sequence in wt_sequence:
                continue
            _, frameshift_status = key.rsplit('.', 1)
            if frameshift_status == 'inframe_splice_site':
                final_mt_sequence, final_wt_sequence = get_mutated_peptide_with_flanking_sequence(wt_sequence, mt_sequence, min(self.flanking_sequence_length, len(wt_sequence)-1, len(mt_sequence)-1))
            elif frameshift_status == 'frameshift_splice_site':
                final_mt_sequence, final_wt_sequence = get_mutated_frameshift_peptide_with_flanking_sequence(wt_sequence, mt_sequence, min(self.flanking_sequence_length, len(wt_sequence)-1, len(mt_sequence)-1))
            else:
                raise Exception("Unexpected frameshift status {} for record {}. Skipping".format(frameshift_status, identifier))
            if final_mt_sequence and final_wt_sequence:
                output_records.append(SeqRecord(final_mt_sequence, id=f"ALT.{key}", description=""))
                if not self.mutant_only:
                    output_records.append(SeqRecord(final_wt_sequence, id=f"WT.{key}", description=""))

        SeqIO.write(output_records, self.trimmed_fasta_file_path, "fasta")

class PvacfuseGenerateProteinFasta(GenerateProteinFasta):
    def __init__(self, **kwargs):
        self.input = kwargs.pop('input', None)
        self.ref_fasta = kwargs.pop('ref_fasta', None)
        self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', 1000)
        super().__init__(**kwargs)

    def generate_fasta(self):
        from pvactools.lib.fusion_to_kmer_pipeline import FusionToKmerPipeline
        params = {
            'input_file': self.input,
            'output_dir': self.temp_dir,
            'transcript_fasta': self.ref_fasta
        }
        pipeline = FusionToKmerPipeline(**params)
        pipeline.generate_fasta()

    def trim_sequences(self):
        print("Trimming Variant Peptide FASTA")
        records = {}
        keys = set()
        for record in SeqIO.parse(self.fasta_file_path, "fasta"):
            records[record.id] = str(record.seq)
            keys.add(record.id.split('.', 1)[1])

        output_records = []
        for key in sorted(keys, key=lambda x: int(x.split('.', 1)[0])):
            mt_seq = records[f"MT.{key}"]
            wt5_seq = records[f"WT5.{key}"]
            if f"WT3.{key}" in records:
                wt3_seq = records[f"WT3.{key}"]
            else:
                wt3_seq = None
            _, variant_type, position = key.rsplit('.', 2)
            position = int(position)
            start_position = position - self.flanking_sequence_length
            if start_position < 0:
                start_position = 0
            if variant_type == 'frameshift_fusion':
                if self.downstream_sequence_length is None:
                    trimmed_mt_seq = mt_seq[start_position:]
                else:
                    trimmed_mt_seq = mt_seq[start_position:(position + self.downstream_sequence_length)]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    trimmed_wt5_seq = wt5_seq[start_position:(position + self.flanking_sequence_length)]
                    output_records.append(SeqRecord(Seq(trimmed_wt5_seq), id=f"WT5.{key}", description=""))
            else:
                end_position = position + self.flanking_sequence_length
                trimmed_mt_seq = mt_seq[start_position:end_position]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    trimmed_wt5_seq = wt5_seq[start_position:end_position]
                    output_records.append(SeqRecord(Seq(trimmed_wt5_seq), id=f"WT5.{key}", description=""))
                    wt3_position = len(wt3_seq) - len(mt_seq[position:])
                    wt3_start_position = wt3_position - self.flanking_sequence_length
                    if wt3_start_position < 0:
                        wt3_start_position = 0
                    trimmed_wt3_seq = wt3_seq[wt3_start_position:(wt3_position + self.flanking_sequence_length)]
                    output_records.append(SeqRecord(Seq(trimmed_wt3_seq), id=f"WT3.{key}", description=""))

        SeqIO.write(output_records, self.trimmed_fasta_file_path, "fasta")
        print("Completed")
