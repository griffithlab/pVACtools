from pvactools.lib.run_pipeline import RunPipeline
from pvactools.lib.fusion_to_kmer_pipeline import FusionToKmerPipeline
from pvactools.lib.pvacfuse_prediction_pipeline import PvacfusePredictionPipeline
from pvactools.tools.pvacfuse.generate_protein_fasta import PvacfuseGenerateProteinFasta
from pvactools.lib.post_processor import PostProcessor

class PvacfuseRunPipeline(RunPipeline):
    def check_tumor_purity_argument(self):
        pass

    def call_input_to_kmer_pipeline(self):
        params = {
            'output_dir'              : self.base_output_dir,
            'input_file'              : self.input_file,
            'sample_name'             : self.sample_name,
            'transcript_fasta'        : self.ref_fasta,
            'starfusion_file'         : self.starfusion_file,
            'class_i_epitope_length'  : self.class_i_epitope_length,
            'class_ii_epitope_length' : self.class_ii_epitope_length,
            'class_i_hla'             : self.class_i_alleles,
            'class_ii_hla'            : self.class_ii_alleles,
        }
        input_to_kmer_pipeline = FusionToKmerPipeline(**params)
        input_to_kmer_pipeline.execute()
        self.transcript_fasta = input_to_kmer_pipeline.create_file_path('fasta')

    def run_prediction_pipeline(self, params):
        self.predictor = PvacfusePredictionPipeline(**params)
        self.predictor.execute()

    def call_generate_protein_fasta(self, params):
        PvacfuseGenerateProteinFasta(**params).trim_sequences()

    def call_post_processor(self, all_epitopes_file, filtered_file, post_processing_params):
        post_processing_params['input_file'] = all_epitopes_file
        post_processing_params['filtered_report_file'] = filtered_file
        post_processing_params['run_coverage_filter'] = True
        post_processing_params['run_transcript_support_level_filter'] = False
        post_processing_params['run_manufacturability_metrics'] = True
        post_processing_params['run_net_chop'] = True if post_processing_params['net_chop_method'] else False
        post_processing_params['run_netmhc_stab'] = True if post_processing_params['netmhc_stab'] else False
        post_processing_params['file_type'] = 'pVACfuse'
        PostProcessor(**post_processing_params).execute()

    def call_ml_predictor(self):
        pass
