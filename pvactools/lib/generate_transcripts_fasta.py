import tempfile
import os
import shutil
import csv
import re
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class GenerateTranscriptsFasta:
    def __init__(self, **kwargs):
        self.sample_name = kwargs.pop('sample_name', 'tmp')
        if self.sample_name is None:
            self.sample_name = 'tmp'
        self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', 1000)
        self.pass_only = kwargs.pop('pass_only', False)
        self.biotypes = kwargs.pop('biotypes', ['protein_coding'])
        self.allow_incomplete_transcripts = kwargs.pop('allow_incomplete_transcripts', False)
        self.input_tsv = kwargs.pop('input_tsv', None)
        self.output_file = kwargs.pop('output_file', None)
        self.temp_dir = tempfile.mkdtemp()
        self.fasta_file_path = kwargs.pop('fasta_file_path', os.path.join(self.temp_dir, f"{self.sample_name}.transcripts.fa"))

    def execute(self):
        self.generate_fasta()
        shutil.copy(self.fasta_file_path, self.output_file)
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def generate_fasta(self):
        raise Exception("Implement in child class")

class PvacseqGenerateTranscriptsFasta(GenerateTranscriptsFasta):
    def __init__(self, **kwargs):
        self.input_vcf = kwargs.pop('input_vcf', None)
        self.phased_proximal_variants_vcf = kwargs.pop('phased_proximal_variants_vcf', None)
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
        }
        pipeline = VariantToKmerPipeline(**params)
        pipeline.generate_fasta()

class PvacspliceGenerateTranscriptsFasta(GenerateTranscriptsFasta):
    def __init__(self, **kwargs):
        self.input_file = kwargs.pop('input_file', None)
        self.annotated_vcf = kwargs.pop('annotated_vcf', None)
        self.ref_fasta = kwargs.pop('ref_fasta', None)
        self.gtf_file = kwargs.pop('gtf_file', None)
        self.junction_score = kwargs.pop('junction_score', 10)
        self.variant_distance = kwargs.pop('variant_distance', 100)
        self.anchor_types = kwargs.pop('anchor_types', ['A', 'D', 'NDA'])
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
            'downstream_sequence_length'       : self.downstream_sequence_length,
            'normal_sample_name'               : None,
            'keep_tmp_files'                   : False,
            'class_i_epitope_length'           : [],
            'class_ii_epitope_length'          : [],
            'class_i_hla'                      : [],
            'class_ii_hla'                     : [],
        }

        pipeline = JunctionToKmerPipeline(**junction_arguments)
        pipeline.generate_fasta()

class PvacfuseGenerateTranscriptsFasta(GenerateTranscriptsFasta):
    def __init__(self, **kwargs):
        self.input = kwargs.pop('input', None)
        self.ref_fasta = kwargs.pop('ref_fasta', None)
        super().__init__(**kwargs)

    def generate_fasta(self):
        from pvactools.lib.fusion_to_kmer_pipeline import FusionToKmerPipeline
        params = {
            'input_file': self.input,
            'output_dir': self.temp_dir,
            'sample_name': self.sample_name,
            'transcript_fasta': self.ref_fasta,
            'downstream_sequence_length': self.downstream_sequence_length,
        }
        pipeline = FusionToKmerPipeline(**params)
        pipeline.generate_fasta()
