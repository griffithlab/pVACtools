import os

from pvactools.lib.run_pipeline import RunPipeline
from pvactools.lib.sequence_fasta_to_vector_fasta import SequenceFastaToVectorFasta
from pvactools.lib.fasta_to_kmers import SequenceFastaToKmers
from pvactools.lib.pvacbind_prediction_pipeline import PvacbindPredictionPipeline

class PvacvectorRunPipeline(RunPipeline):
    def check_tumor_purity_argument(self):
        pass

    def check_downstream_sequence_length_argument(self):
        pass

    def call_input_to_kmer_pipeline(self):
        for length in set(self.class_i_epitope_length + self.class_ii_epitope_length):
            tmp_dir = os.path.join(self.output_dir, 'tmp')
            os.makedirs(tmp_dir, exist_ok=True)
            output_file = os.path.join(tmp_dir, f"{self.sample_name}.{length}.fa")
            params = {
                'input_file': self.input_file,
                'output_file': output_file,
                'epitope_length': length,
                'spacer': self.spacer,
                'junctions_to_test': self.junctions_to_test,
                'clip_length': self.clip_length,
            }
            fasta_generator = SequenceFastaToVectorFasta(**params)
            fasta_generator.execute()

            fasta_to_kmer_arguments = {
                'fasta': output_file,
                'output_dir': self.output_dir,
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
        pass

    def call_ml_predictor(self):
        pass
