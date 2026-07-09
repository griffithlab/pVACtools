from pvactools.lib.run_pipeline import RunPipeline
from pvactools.lib.fasta_to_kmers import SequenceFastaToKmers
from pvactools.lib.pvacbind_prediction_pipeline import PvacbindPredictionPipeline
from pvactools.lib.post_processor import PostProcessor

class PvacbindRunPipeline(RunPipeline):
    def check_tumor_purity_argument(self):
        pass

    def check_downstream_sequence_length_argument(self):
        pass

    def call_input_to_kmer_pipeline(self):
        for length in set(self.class_i_epitope_length + self.class_ii_epitope_length):
            fasta_to_kmer_arguments = {
                'fasta': self.input_file,
                'output_dir': self.base_output_dir,
                'epitope_length': length,
                'sample_name': self.sample_name,
            }
            SequenceFastaToKmers(**fasta_to_kmer_arguments).execute()
        self.transcript_fasta = self.input_file

    def run_prediction_pipeline(self, params):
        self.predictor = PvacbindPredictionPipeline(**params)
        self.predictor.execute()

    def generate_flanked_protein_fasta(self, flanking_sequence_length):
        return self.input_file

    def call_generate_protein_fasta(self, params):
        pass

    def call_post_processor(self, all_epitopes_file, filtered_file, post_processing_params):
        post_processing_params['input_file'] = all_epitopes_file
        post_processing_params['filtered_report_file'] = filtered_file
        post_processing_params['run_coverage_filter'] = False
        post_processing_params['run_transcript_support_level_filter'] = False
        post_processing_params['run_manufacturability_metrics'] = True
        post_processing_params['run_net_chop'] = True if post_processing_params['net_chop_method'] else False
        post_processing_params['run_netmhc_stab'] = True if post_processing_params['netmhc_stab'] else False
        post_processing_params['file_type'] = 'pVACbind'
        PostProcessor(**post_processing_params).execute()

    def call_ml_predictor(self):
        pass
