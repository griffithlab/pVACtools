import os
import logging

from pvactools.lib.run_pipeline import RunPipeline
from pvactools.lib.variant_to_kmer_pipeline import VariantToKmerPipeline
from pvactools.lib.pvacseq_prediction_pipeline import PvacseqPredictionPipeline
from pvactools.tools.pvacseq.generate_protein_fasta import PvacseqGenerateProteinFasta
from pvactools.lib.post_processor import PvacseqPostProcessor

class PvacseqRunPipeline(RunPipeline):
    def call_input_to_kmer_pipeline(self):
        params = {
            'output_dir'                  : self.base_output_dir,
            'input_file'                  : self.input_file,
            'sample_name'                 : self.sample_name,
            'pass_only'                   : self.pass_only,
            'normal_sample_name'          : self.normal_sample_name,
            'proximal_variants_vcf'       : self.phased_proximal_variants_vcf,
            'biotypes'                    : self.biotypes,
            'allow_incomplete_transcripts': self.allow_incomplete_transcripts,
            'downstream_sequence_length'  : self.downstream_sequence_length,
            'class_i_epitope_length'      : self.class_i_epitope_length,
            'class_ii_epitope_length'     : self.class_ii_epitope_length,
            'class_i_hla'                 : self.class_i_alleles,
            'class_ii_hla'                : self.class_ii_alleles,
        }
        input_to_kmer_pipeline = VariantToKmerPipeline(**params)
        input_to_kmer_pipeline.execute()
        self.transcript_fasta = input_to_kmer_pipeline.create_file_path('fasta')

    def run_prediction_pipeline(self, params):
        self.predictor = PvacseqPredictionPipeline(**params)
        self.predictor.execute()

    def call_generate_protein_fasta(self, params):
        PvacseqGenerateProteinFasta(**params).trim_sequences()

    def call_post_processor(self, all_epitopes_file, filtered_file, post_processing_params):
        post_processing_params['input_file'] = all_epitopes_file
        post_processing_params['filtered_report_file'] = filtered_file
        PvacseqPostProcessor(**post_processing_params).execute()

    def call_ml_predictor(self):
        if len(self.class_i_prediction_algorithms) > 0 and len(self.class_i_alleles) > 0 and len(self.class_ii_prediction_algorithms) > 0 and len(self.class_ii_alleles) > 0:
            # Run ML predictions
            if self.run_ml_predictions:
                logging.info("Running ML predictions...")
                if not 'all' in self.prediction_algorithms:
                    logging.info("Caution: Use 'all' in prediction_algorithms is strongly recommended. Missing features will be filled with NA and will cause predictions to be inaccurate. Running ML predictions regardless...")

                # Locate input files
                file1 = os.path.join(self.base_output_dir, 'MHC_Class_I', f"{self.sample_name}.MHC_I.all_epitopes.aggregated.tsv")
                file2 = os.path.join(self.base_output_dir, 'MHC_Class_I', f"{self.sample_name}.MHC_I.all_epitopes.tsv")
                file3 = os.path.join(self.base_output_dir, 'MHC_Class_II', f"{self.sample_name}.MHC_II.all_epitopes.aggregated.tsv")
                file4 = os.path.join(self.base_output_dir, 'MHC_Class_I', f"{self.sample_name}.MHC_I.all_epitopes.aggregated.metrics.json")

                # Check if all required files exist
                required_files = [file1, file2, file3, file4]
                missing_files = [f for f in required_files if not os.path.exists(f)]
                if missing_files:
                    logging.warning(f"Warning: Missing required files for ML predictions: {missing_files}")
                    logging.warning("Skipping ML predictions.")
                    return

                # Save ML output in the same folder as MHC_I.all_epitopes.aggregated.tsv (MHC_Class_I)
                ml_output_dir = os.path.dirname(file1)

                try:
                    # Import and run ML predictions
                    from pvactools.lib.ml_predictor import run_ml_predictions

                    output_file = run_ml_predictions(
                        class1_aggregated_path=file1,
                        class1_all_epitopes_path=file2,
                        class2_aggregated_path=file3,
                        model_artifacts_path=None,  # None uses default package location
                        output_dir=ml_output_dir,
                        sample_name=self.sample_name,
                        ml_threshold_accept=self.ml_threshold_accept,
                        ml_threshold_reject=self.ml_threshold_reject
                    )
                    logging.info(f"ML predictions completed successfully using Class I and Class II files. Results saved to: {output_file}")

                except Exception as e:
                    logging.warning(f"Error during standalone ML predictions: {str(e)}")
                    logging.warning("Continuing with pipeline without ML predictions.")
