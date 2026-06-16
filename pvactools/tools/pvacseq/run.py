import sys
import argparse
import os
import shutil
import yaml
import platform
import copy

from pvactools.lib.prediction_class import *
from pvactools.lib.pipeline import PvacsplicePipeline
from pvactools.lib.run_argument_parser import PvacseqRunArgumentParser
from pvactools.lib.post_processor import PostProcessor
from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *
from pvactools.lib.variant_pipeline import VariantPipeline
from pvactools.tools.pvacseq.generate_protein_fasta import PvacseqGenerateProteinFasta
from pvactools.lib.print_log import *

def define_parser():
    return PvacseqRunArgumentParser().parser

def create_per_class_report(files, all_epitopes_output_file, filtered_report_file, post_processing_params, run_params):
    for file_name in files:
        if not os.path.exists(file_name):
            print("File {} doesn't exist. Aborting.".format(file_name))
            return

    combine_reports(files, all_epitopes_output_file)

    post_processing_params['input_file'] = all_epitopes_output_file
    post_processing_params['filtered_report_file'] = filtered_report_file
    post_processing_params['minimum_fold_change'] = None
    post_processing_params['run_coverage_filter'] = True
    post_processing_params['run_transcript_support_level_filter'] = False
    post_processing_params['run_manufacturability_metrics'] = True
    if run_params['net_chop_method']:
        post_processing_params['run_net_chop'] = True
        post_processing_params['net_chop_fasta'] = run_params['net_chop_fasta']
    else:
        post_processing_params['run_net_chop'] = False
    post_processing_params['run_netmhc_stab'] = True if run_params['netmhc_stab'] else False
    if run_params['run_reference_proteome_similarity']:
        post_processing_params['fasta'] = run_params['fasta']
    post_processing_params['species'] = run_params['species']
    post_processing_params['file_type'] = 'pVACseq'
    PostProcessor(**post_processing_params).execute()

def create_combined_reports(base_output_dir, args):
    output_dir = os.path.join(base_output_dir, 'combined')
    os.makedirs(output_dir, exist_ok=True)

    file1 = os.path.join(base_output_dir, 'MHC_Class_I', "{}.MHC_I.all_epitopes.tsv".format(args.sample_name))
    file2 = os.path.join(base_output_dir, 'MHC_Class_II', "{}.MHC_II.all_epitopes.tsv".format(args.sample_name))
    if not os.path.exists(file1):
        print("File {} doesn't exist. Aborting.".format(file1))
        return
    if not os.path.exists(file2):
        print("File {} doesn't exist. Aborting.".format(file2))
        return

    combined_output_file = os.path.join(output_dir, "{}.Combined.all_epitopes.tsv".format(args.sample_name))
    combine_reports([file1, file2], combined_output_file)
    filtered_report_file = os.path.join(output_dir, "{}.Combined.filtered.tsv".format(args.sample_name))

    post_processing_params = vars(args)
    post_processing_params['input_file'] = combined_output_file
    post_processing_params['filtered_report_file'] = filtered_report_file
    post_processing_params['run_coverage_filter'] = True
    post_processing_params['run_transcript_support_level_filter'] = True
    post_processing_params['run_net_chop'] = False
    post_processing_params['run_netmhc_stab'] = False
    post_processing_params['run_manufacturability_metrics'] = False
    post_processing_params['run_reference_proteome_similarity'] = False
    post_processing_params['file_type'] = 'pVACseq'
    post_processing_params['filename_addition'] = "Combined"

    PostProcessor(**post_processing_params).execute()

def locate_ml_input_files(base_output_dir, sample_name):
    """
    Locate the three input files required for ML predictions.

    Args:
        base_output_dir (str): Base output directory
        sample_name (str): Sample name

    Returns:
        tuple: Paths to the three required files (file1, file2, file3)
    """
    file1 = os.path.join(base_output_dir, 'MHC_Class_I', "{}.MHC_I.all_epitopes.aggregated.tsv".format(sample_name))
    file2 = os.path.join(base_output_dir, 'MHC_Class_I', "{}.MHC_I.all_epitopes.tsv".format(sample_name))
    file3 = os.path.join(base_output_dir, 'MHC_Class_II', "{}.MHC_II.all_epitopes.aggregated.tsv".format(sample_name))
    file4 = os.path.join(base_output_dir, 'MHC_Class_I', "{}.MHC_I.all_epitopes.aggregated.metrics.json".format(sample_name))

    return file1, file2, file3, file4

def run_ml_predictions(base_output_dir, args):
    """
    Run ML predictions as a standalone process when both Class I and Class II predictions are available.

    Args:
        base_output_dir (str): Base output directory
        args: Command line arguments
    """

    print("Running ML predictions...")
    if not 'all' in args.prediction_algorithms:
        print("Caution: Use 'all' in prediction_algorithms is strongly recommended. Missing features will be filled with NA and will cause predictions to be inaccurate. Running ML predictions regardless...")

    # Locate input files
    file1, file2, file3, file4 = locate_ml_input_files(base_output_dir, args.sample_name)

    # Check if all required files exist
    required_files = [file1, file2, file3, file4]
    missing_files = [f for f in required_files if not os.path.exists(f)]
    if missing_files:
        print(f"Warning: Missing required files for ML predictions: {missing_files}")
        print("Skipping ML predictions.")
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
            sample_name=args.sample_name,
            ml_threshold_accept=args.ml_threshold_accept,
            ml_threshold_reject=args.ml_threshold_reject
        )
        print(f"ML predictions completed successfully using Class I and Class II files. Results saved to: {output_file}")

    except Exception as e:
        print(f"Error during standalone ML predictions: {str(e)}")
        print("Continuing with pipeline without ML predictions.")

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.fasta_size%2 != 0:
        sys.exit("The fasta size needs to be an even number")

    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")

    if args.downstream_sequence_length == 'full':
        downstream_sequence_length = None
    elif args.downstream_sequence_length.isdigit():
        downstream_sequence_length = int(args.downstream_sequence_length)
    else:
        sys.exit("The downstream sequence length needs to be a positive integer or 'full'")

    if args.tumor_purity is not None:
        if args.tumor_purity > 1:
            raise Exception("--tumor-purity must be a float between 0 and 1. Value too large: {}".format(args.tumor_purity))
        elif args.tumor_purity < 0:
            raise Exception("--tumor-purity must be a float between 0 and 1. Value too small: {}".format(args.tumor_purity))

    if args.n_threads > 1 and platform.system() == "Darwin":
        raise Exception("Multithreading is not supported on MacOS")

    if (args.netmhciipan_version == '4.0' and args.iedb_install_directory is not None):
        raise Exception("Standalone IEDB does not support version 4.0")
    NetMHCIIVersion.netmhciipan_version = args.netmhciipan_version

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    alleles = combine_class_ii_alleles(args.allele)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(alleles)

    input_file_type = 'vcf'
    base_output_dir = os.path.abspath(args.output_dir)
    os.makedirs(base_output_dir, exist_ok=True)

    print_log(os.path.join(base_output_dir, 'log'), vars(args), 'inputs')

    variant_arguments = {
        'output_dir'                  : base_output_dir,
        'input_file'                  : args.input_file,
        'sample_name'                 : args.sample_name,
        'pass_only'                   : args.pass_only,
        'normal_sample_name'          : args.normal_sample_name,
        'proximal_variants_vcf'       : args.phased_proximal_variants_vcf,
        'biotypes'                    : args.biotypes,
        'allow_incomplete_transcripts': args.allow_incomplete_transcripts,
        'downstream_sequence_length'  : downstream_sequence_length,
        'class_i_epitope_length'      : args.class_i_epitope_length,
        'class_ii_epitope_length'     : args.class_ii_epitope_length,
        'class_i_hla'                 : class_i_alleles,
        'class_ii_hla'                : class_ii_alleles,
    }
    variant_pipeline = VariantPipeline(**variant_arguments)
    variant_pipeline.execute()

    shared_arguments = {
        'input_file'                : args.input_file,
        'input_file_type'           : input_file_type,
        'sample_name'               : args.sample_name,
        'top_score_metric'          : args.top_score_metric,
        'top_score_metric2'         : args.top_score_metric2,
        'binding_threshold'         : args.binding_threshold,
        'binding_percentile_threshold': args.binding_percentile_threshold,
        'immunogenicity_percentile_threshold': args.immunogenicity_percentile_threshold,
        'presentation_percentile_threshold': args.presentation_percentile_threshold,
        'percentile_threshold_strategy': args.percentile_threshold_strategy,
        'allele_specific_binding_thresholds': args.allele_specific_binding_thresholds,
        'minimum_fold_change'       : args.minimum_fold_change,
        'net_chop_method'           : args.net_chop_method,
        'net_chop_threshold'        : args.net_chop_threshold,
        'normal_cov'                : args.normal_cov,
        'normal_vaf'                : args.normal_vaf,
        'tdna_cov'                  : args.tdna_cov,
        'tdna_vaf'                  : args.tdna_vaf,
        'trna_cov'                  : args.trna_cov,
        'trna_vaf'                  : args.trna_vaf,
        'expn_val'                  : args.expn_val,
        'additional_report_columns' : args.additional_report_columns,
        'fasta_size'                : args.fasta_size,
        'iedb_retries'              : args.iedb_retries,
        'downstream_sequence_length': downstream_sequence_length,
        'keep_tmp_files'            : args.keep_tmp_files,
        'pass_only'                 : args.pass_only,
        'normal_sample_name'        : args.normal_sample_name,
        'phased_proximal_variants_vcf' : args.phased_proximal_variants_vcf,
        'n_threads'                 : args.n_threads,
        'transcript_prioritization_strategy': args.transcript_prioritization_strategy,
        'maximum_transcript_support_level': args.maximum_transcript_support_level,
        'biotypes'                  : args.biotypes,
        'allow_incomplete_transcripts': args.allow_incomplete_transcripts,
        'species'                   : species,
        'run_reference_proteome_similarity': args.run_reference_proteome_similarity,
        'blastp_path'               : args.blastp_path,
        'blastp_db'                 : args.blastp_db,
        'tumor_purity'              : args.tumor_purity,
        'problematic_amino_acids'   : args.problematic_amino_acids,
        'peptide_fasta'             : args.peptide_fasta,
        'allele_specific_anchors'   : args.allele_specific_anchors,
        'anchor_contribution_threshold' : args.anchor_contribution_threshold,
        'aggregate_inclusion_binding_threshold': args.aggregate_inclusion_binding_threshold,
        'aggregate_inclusion_count_limit': args.aggregate_inclusion_count_limit,
        'genes_of_interest_file': args.genes_of_interest_file,
        'run_ml_predictions': args.run_ml_predictions,
        'ml_threshold_accept': args.ml_threshold_accept,
        'ml_threshold_reject': args.ml_threshold_reject,
    }

    if args.iedb_install_directory:
        iedb_mhc_i_executable = os.path.join(args.iedb_install_directory, 'mhc_i', 'src', 'predict_binding.py')
        if not os.path.exists(iedb_mhc_i_executable):
            sys.exit("IEDB MHC I executable path doesn't exist %s" % iedb_mhc_i_executable)
        iedb_mhc_ii_executable = os.path.join(args.iedb_install_directory, 'mhc_ii', 'mhc_II_binding.py')
        if not os.path.exists(iedb_mhc_ii_executable):
            sys.exit("IEDB MHC II executable path doesn't exist %s" % iedb_mhc_ii_executable)
    else:
        iedb_mhc_i_executable = None
        iedb_mhc_ii_executable = None

    if args.use_normalized_percentiles and species != 'human':
        print("WARNING: Normalized percentiles are only available for human alleles. Option will be ignored.")
        args.use_normalized_percentiles = False

    all_params = {
        'I': {
            'iedb_executable': iedb_mhc_i_executable,
            'prediction_algorithms': class_i_prediction_algorithms,
            'alleles': class_i_alleles,
            'epitope_lengths': args.class_i_epitope_length,
            'netmhc_stab': args.netmhc_stab,
            'use_normalized_percentiles': args.use_normalized_percentiles,
            'reference_scores_path': args.reference_scores_path
        },
        'II': {
            'iedb_executable': iedb_mhc_ii_executable,
            'prediction_algorithms': class_ii_prediction_algorithms,
            'alleles': class_ii_alleles,
            'epitope_lengths': args.class_ii_epitope_length,
            'netmhc_stab': False,
            'use_normalized_percentiles': False,
            'reference_scores_path': args.reference_scores_path
        }
    }

    for (mhc_class, params) in all_params.items():
        prediction_algorithms = params['prediction_algorithms']
        alleles = params['alleles']
        epitope_lengths = params['epitope_lengths']
        iedb_executable = params['iedb_executable']
        netmhc_stab = params['netmhc_stab']
        use_normalized_percentiles = params['use_normalized_percentiles']
        reference_scores_path = params['reference_scores_path']

        if len(prediction_algorithms) > 0 and len(alleles) > 0:
            print("Executing MHC Class {} predictions".format(mhc_class))

            output_dir = os.path.join(base_output_dir, 'MHC_Class_{}'.format(mhc_class))
            os.makedirs(output_dir, exist_ok=True)

            output_files = []
            run_arguments = copy.deepcopy(shared_arguments)
            run_arguments['alleles']               = alleles
            run_arguments['iedb_executable']       = iedb_executable
            run_arguments['prediction_algorithms'] = prediction_algorithms
            run_arguments['netmhc_stab']           = netmhc_stab
            run_arguments['use_normalized_percentiles'] = use_normalized_percentiles
            run_arguments['reference_scores_path'] = reference_scores_path

            for epitope_length in epitope_lengths:
                per_length_run_arguments = copy.deepcopy(run_arguments)
                per_epitope_output_dir = os.path.join(output_dir, str(epitope_length))
                os.makedirs(per_epitope_output_dir, exist_ok=True)
                input_file = os.path.join(args.output_dir, f'{args.sample_name}.{epitope_length}.fa')
                if os.path.getsize(input_file) == 0:
                    print("The intermediate FASTA file for epitope length {} is empty. No processable variants found.")
                    continue

                per_length_run_arguments['input_file']      = input_file
                per_length_run_arguments['epitope_lengths'] = epitope_length
                per_length_run_arguments['output_dir']      = per_epitope_output_dir
                pipeline = PvacsplicePipeline(**per_length_run_arguments)
                pipeline.execute()
                output_file = os.path.join(per_epitope_output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
                if os.path.exists(output_file):
                    output_files.append(output_file)
            if len(output_files) > 0:
                # copy fasta to output dir
                transcript_fasta = variant_pipeline.create_file_path('fasta')
                fasta_file = os.path.join(output_dir, "{}.fasta".format(args.sample_name))
                params = {
                    'fasta_file_path': transcript_fasta,
                    'trimmed_fasta_file_path': fasta_file,
                    'flanking_sequence_length': max(epitope_lengths) - 1,
                    'mutant_only': False,
                }
                PvacseqGenerateProteinFasta(**params).trim_sequences()
                if args.run_reference_proteome_similarity:
                    fasta_file = os.path.join(output_dir, "{}.7.fasta".format(args.sample_name))
                    if not os.path.exists(fasta_file):
                        params = {
                            'fasta_file_path': transcript_fasta,
                            'trimmed_fasta_file_path': fasta_file,
                            'flanking_sequence_length': 7,
                            'mutant_only': False,
                        }
                        PvacseqGenerateProteinFasta(**params).trim_sequences()
                    run_arguments['fasta'] = fasta_file
                # generate and copy net_chop fasta to output dir if specified
                if args.net_chop_method:
                    fasta_file = os.path.join(output_dir, "{}.10.fasta".format(args.sample_name))
                    if not os.path.exists(fasta_file):
                        params = {
                            'fasta_file_path': transcript_fasta,
                            'trimmed_fasta_file_path': fasta_file,
                            'flanking_sequence_length': 10,
                            'mutant_only': False,
                        }
                        PvacseqGenerateProteinFasta(**params).trim_sequences()
                    run_arguments['net_chop_fasta'] = fasta_file
                all_epitopes_file = os.path.join(output_dir, "{}.MHC_{}.all_epitopes.tsv".format(args.sample_name,mhc_class))
                filtered_file = os.path.join(output_dir, "{}.MHC_{}.filtered.tsv".format(args.sample_name,mhc_class))
                post_processing_params = vars(args).copy()
                post_processing_params["filename_addition"] = "MHC_{}".format(mhc_class)
                create_per_class_report(output_files, all_epitopes_file, filtered_file, post_processing_params, run_arguments)
            else:
                print("\nNo processable variants found. Aborting.\n")
        elif len(prediction_algorithms) == 0:
            print("No MHC class {} prediction algorithms chosen. Skipping MHC class {} predictions.".format(mhc_class, mhc_class))
        elif len(alleles) == 0:
            print("No MHC class {} alleles chosen. Skipping MHC class {} predictions.".format(mhc_class, mhc_class))

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0 and len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        print("Creating combined reports")
        create_combined_reports(base_output_dir, args)

        # Run ML predictions
        if args.run_ml_predictions:
            run_ml_predictions(base_output_dir, args)

    change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
