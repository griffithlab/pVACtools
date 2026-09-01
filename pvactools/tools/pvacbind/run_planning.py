import os
from dataclasses import dataclass


# pVACbind-local run planning, output layout, and handoff construction.
# Execution, parser validation, prediction classification, and report semantics live elsewhere.
MHC_CLASS_I_DIRECTORY = "MHC_Class_I"
MHC_CLASS_II_DIRECTORY = "MHC_Class_II"
COMBINED_DIRECTORY = "combined"

MHC_CLASS_I_FILENAME_ADDITION = "MHC_I"
MHC_CLASS_II_FILENAME_ADDITION = "MHC_II"
COMBINED_FILENAME_ADDITION = "Combined"


@dataclass(frozen=True)
class PvacbindOutputLayout:
    base_output_dir: str
    sample_name: str

    @property
    def class_i_output_dir(self):
        return os.path.join(self.base_output_dir, MHC_CLASS_I_DIRECTORY)

    @property
    def class_ii_output_dir(self):
        return os.path.join(self.base_output_dir, MHC_CLASS_II_DIRECTORY)

    @property
    def combined_output_dir(self):
        return os.path.join(self.base_output_dir, COMBINED_DIRECTORY)

    @property
    def class_i_all_epitopes_file(self):
        return os.path.join(
            self.class_i_output_dir,
            "{}.{}.all_epitopes.tsv".format(self.sample_name, MHC_CLASS_I_FILENAME_ADDITION),
        )

    @property
    def class_ii_all_epitopes_file(self):
        return os.path.join(
            self.class_ii_output_dir,
            "{}.{}.all_epitopes.tsv".format(self.sample_name, MHC_CLASS_II_FILENAME_ADDITION),
        )

    @property
    def combined_all_epitopes_file(self):
        return os.path.join(
            self.combined_output_dir,
            "{}.{}.all_epitopes.tsv".format(self.sample_name, COMBINED_FILENAME_ADDITION),
        )

    @property
    def combined_filtered_report_file(self):
        return os.path.join(
            self.combined_output_dir,
            "{}.{}.filtered.tsv".format(self.sample_name, COMBINED_FILENAME_ADDITION),
        )


@dataclass(frozen=True)
class PvacbindRunPlan:
    class_i_prediction_algorithms: list
    class_ii_prediction_algorithms: list
    class_i_alleles: list
    class_ii_alleles: list
    species: str

    def should_run_class_i(self):
        return len(self.class_i_prediction_algorithms) > 0 and len(self.class_i_alleles) > 0

    def should_run_class_ii(self):
        return len(self.class_ii_prediction_algorithms) > 0 and len(self.class_ii_alleles) > 0

    def should_create_combined_reports(self):
        return self.should_run_class_i() and self.should_run_class_ii()

    def class_i_skip_message(self):
        if self.should_run_class_i():
            return None
        elif len(self.class_i_prediction_algorithms) == 0:
            return "No MHC class I prediction algorithms chosen. Skipping MHC class I predictions."
        elif len(self.class_i_alleles) == 0:
            return "No MHC class I alleles chosen. Skipping MHC class I predictions."
        return None

    def class_ii_skip_message(self):
        if self.should_run_class_ii():
            return None
        elif len(self.class_ii_prediction_algorithms) == 0:
            return "No MHC class II prediction algorithms chosen. Skipping MHC class II predictions."
        elif len(self.class_ii_alleles) == 0:
            return "No MHC class II alleles chosen. Skipping MHC class II predictions."
        return None


@dataclass(frozen=True)
class PvacbindPipelineHandoffBuilder:
    args: object
    run_plan: PvacbindRunPlan
    output_layout: PvacbindOutputLayout

    def shared_arguments(self):
        return {
            'input_file'                : self.args.input_file,
            'input_file_type'           : 'fasta',
            'sample_name'               : self.args.sample_name,
            'top_score_metric'          : self.args.top_score_metric,
            'top_score_metric2'         : self.args.top_score_metric2,
            'binding_threshold'         : self.args.binding_threshold,
            'binding_percentile_threshold': self.args.binding_percentile_threshold,
            'immunogenicity_percentile_threshold': self.args.immunogenicity_percentile_threshold,
            'presentation_percentile_threshold': self.args.presentation_percentile_threshold,
            'percentile_threshold_strategy': self.args.percentile_threshold_strategy,
            'allele_specific_binding_thresholds': self.args.allele_specific_binding_thresholds,
            'net_chop_fasta'            : self.args.input_file,
            'net_chop_method'           : self.args.net_chop_method,
            'net_chop_threshold'        : self.args.net_chop_threshold,
            'additional_report_columns' : self.args.additional_report_columns,
            'fasta_size'                : self.args.fasta_size,
            'iedb_retries'              : self.args.iedb_retries,
            'keep_tmp_files'            : self.args.keep_tmp_files,
            'n_threads'                 : self.args.n_threads,
            'species'                   : self.run_plan.species,
            'run_reference_proteome_similarity': self.args.run_reference_proteome_similarity,
            'blastp_path'               : self.args.blastp_path,
            'blastp_db'                 : self.args.blastp_db,
            'problematic_amino_acids'   : self.args.problematic_amino_acids,
            'run_post_processor'        : True,
            'peptide_fasta'             : self.args.peptide_fasta,
            'aggregate_inclusion_binding_threshold': self.args.aggregate_inclusion_binding_threshold,
            'aggregate_inclusion_count_limit': self.args.aggregate_inclusion_count_limit,
        }

    def class_i_arguments(self, iedb_mhc_i_executable):
        class_i_arguments = self.shared_arguments()
        class_i_arguments['alleles']                 = self.run_plan.class_i_alleles
        class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
        class_i_arguments['epitope_lengths']         = self.args.class_i_epitope_length
        class_i_arguments['prediction_algorithms']   = self.run_plan.class_i_prediction_algorithms
        class_i_arguments['output_dir']              = self.output_layout.class_i_output_dir
        class_i_arguments['netmhc_stab']             = self.args.netmhc_stab
        class_i_arguments['filename_addition']       = MHC_CLASS_I_FILENAME_ADDITION
        class_i_arguments['use_normalized_percentiles'] = self.args.use_normalized_percentiles
        class_i_arguments['reference_scores_path']   = self.args.reference_scores_path
        return class_i_arguments

    def class_ii_arguments(self, iedb_mhc_ii_executable):
        class_ii_arguments = self.shared_arguments()
        class_ii_arguments['alleles']                 = self.run_plan.class_ii_alleles
        class_ii_arguments['prediction_algorithms']   = self.run_plan.class_ii_prediction_algorithms
        class_ii_arguments['iedb_executable']         = iedb_mhc_ii_executable
        class_ii_arguments['epitope_lengths']         = self.args.class_ii_epitope_length
        class_ii_arguments['output_dir']              = self.output_layout.class_ii_output_dir
        class_ii_arguments['netmhc_stab']             = False
        class_ii_arguments['filename_addition']       = MHC_CLASS_II_FILENAME_ADDITION
        return class_ii_arguments


@dataclass(frozen=True)
class CombinedReportHandoff:
    output_layout: PvacbindOutputLayout

    @property
    def input_files(self):
        return [
            self.output_layout.class_i_all_epitopes_file,
            self.output_layout.class_ii_all_epitopes_file,
        ]

    def post_processing_params(self, args):
        post_processing_params = vars(args).copy()
        post_processing_params['input_file'] = self.output_layout.combined_all_epitopes_file
        post_processing_params['filtered_report_file'] = self.output_layout.combined_filtered_report_file
        post_processing_params['run_coverage_filter'] = False
        post_processing_params['minimum_fold_change'] = None
        post_processing_params['file_type'] = 'pVACbind'
        post_processing_params['run_transcript_support_level_filter'] = False
        post_processing_params['run_net_chop'] = False
        post_processing_params['run_netmhc_stab'] = False
        post_processing_params['run_manufacturability_metrics'] = False
        post_processing_params['run_reference_proteome_similarity'] = False
        post_processing_params["filename_addition"] = COMBINED_FILENAME_ADDITION
        return post_processing_params
