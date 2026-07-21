from pvactools.lib.run_pipeline import RunPipeline
from pvactools.lib.fasta_to_kmers import SequenceFastaToKmers
from pvactools.lib.pvacbind_prediction_pipeline import PvacbindPredictionPipeline
from pvactools.lib.post_processor import PvacbindPostProcessor

class PvacbindRunPipeline(RunPipeline):
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
        PvacbindPostProcessor(**post_processing_params).execute()

    def call_ml_predictor(self):
        pass
