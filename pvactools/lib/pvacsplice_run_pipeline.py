from pathlib import Path
import shutil

from pvactools.lib.run_utils import *
from pvactools.lib.run_pipeline import RunPipeline
from pvactools.lib.junction_to_kmer_pipeline import JunctionToKmerPipeline
from pvactools.lib.pvacsplice_prediction_pipeline import PvacsplicePredictionPipeline
from pvactools.lib.generate_protein_fasta import PvacspliceGenerateProteinFasta
from pvactools.lib.post_processor import PvacsplicePostProcessor

class PvacspliceRunPipeline(RunPipeline):
    def check_downstream_sequence_length_argument(self):
        pass

    def extra_argument_checks(self):
        # ref fasta
        if Path(self.ref_fasta).suffix not in ['.fa', '.fasta']:
            raise Exception('The fasta input path does not point to a fasta file.')
        if is_gz_file(Path(self.ref_fasta)):
            raise Exception('pVACsplice does not currently support gzipped fasta files.')

        # gtf
        if Path(self.gtf_file).suffix not in ['.gtf', '.tsv'] and not is_gz_file(self.gtf_file):
            raise Exception('The gtf input path does not point to a gtf file.')

        # vcf
        if Path(self.annotated_vcf).suffix != '.vcf' and not is_gz_file(self.annotated_vcf):
            raise Exception('The vcf input path does not point to a vcf file.')

        # vcf gz.tbi index file
        if is_gz_file(self.annotated_vcf) and not Path(f'{self.annotated_vcf}.tbi').exists():
            raise Exception('Gzipped VCF files must be indexed. (tabix -p vcf <vcf_file>)')

    def call_input_to_kmer_pipeline(self):
        params = {
            'input_file_type'              : 'junctions',
            'junctions_dir'                : self.base_output_dir,
            'input_file'                   : self.input_file,
            'gtf_file'                     : self.gtf_file,
            'save_gtf'                     : self.save_gtf,
            'sample_name'                  : self.sample_name,
            'ref_fasta'                    : self.ref_fasta,
            'annotated_vcf'                : self.annotated_vcf,
            'pass_only'                    : self.pass_only,
            'class_i_epitope_length'       : self.class_i_epitope_length,
            'class_ii_epitope_length'      : self.class_ii_epitope_length,
            'biotypes'                     : self.biotypes,
            'allow_incomplete_transcripts' : self.allow_incomplete_transcripts,
            'junction_score'               : self.junction_score,
            'variant_distance'             : self.variant_distance,
            'anchor_types'                 : self.anchor_types,
            'normal_sample_name'           : self.normal_sample_name,
            'class_i_hla'                  : self.class_i_alleles,
            'class_ii_hla'                 : self.class_ii_alleles,
            'keep_tmp_files'               : self.keep_tmp_files,
        }
        input_to_kmer_pipeline = JunctionToKmerPipeline(**params)
        input_to_kmer_pipeline.execute()
        self.transcript_fasta = input_to_kmer_pipeline.create_file_path('fasta')

    def run_prediction_pipeline(self, params):
        self.predictor = PvacsplicePredictionPipeline(**params)
        self.predictor.execute()

    def call_generate_protein_fasta(self, params):
        PvacspliceGenerateProteinFasta(**params).trim_sequences()

    def call_post_processor(self, all_epitopes_file, filtered_file, post_processing_params):
        post_processing_params['input_file'] = all_epitopes_file
        post_processing_params['filtered_report_file'] = filtered_file
        PvacsplicePostProcessor(**post_processing_params).execute()

    def call_ml_predictor(self):
        pass

    def extra_steps(self):
        if self.save_gtf is False:
            shutil.rmtree(os.path.join(self.base_output_dir, 'tmp'), ignore_errors=True)
