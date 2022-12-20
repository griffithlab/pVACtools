import sys
import argparse
import os
import shutil
import yaml

from pvactools.lib.prediction_class import *
from pvactools.lib.pipeline import PvacbindPipeline
from pvactools.lib.run_argument_parser import PvacbindRunArgumentParser
from pvactools.lib.post_processor import PostProcessor
import pvactools.lib.run_utils

def define_parser():
    return PvacbindRunArgumentParser().parser

def create_combined_reports(base_output_dir, args):
    output_dir = os.path.join(base_output_dir, 'combined')
    os.makedirs(output_dir, exist_ok=True)

    file1 = os.path.join(base_output_dir, 'MHC_Class_I', "{}.all_epitopes.tsv".format(args.sample_name))
    file2 = os.path.join(base_output_dir, 'MHC_Class_II', "{}.all_epitopes.tsv".format(args.sample_name))
    if not os.path.exists(file1):
        print("File {} doesn't exist. Aborting.".format(file1))
        return
    if not os.path.exists(file2):
        print("File {} doesn't exist. Aborting.".format(file2))
        return
    combined_output_file = os.path.join(output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
    pvactools.lib.run_utils.combine_reports([file1, file2], combined_output_file)
    filtered_report_file = os.path.join(output_dir, "{}.filtered.tsv".format(args.sample_name))

    post_processing_params = vars(args)
    post_processing_params['input_file'] = combined_output_file
    post_processing_params['filtered_report_file'] = filtered_report_file
    post_processing_params['run_coverage_filter'] = False
    post_processing_params['minimum_fold_change'] = None
    post_processing_params['file_type'] = 'pVACbind'
    post_processing_params['run_transcript_support_level_filter'] = False
    post_processing_params['run_net_chop'] = False
    post_processing_params['run_netmhc_stab'] = False
    post_processing_params['run_manufacturability_metrics'] = False
    post_processing_params['run_reference_proteome_similarity'] = False

    PostProcessor(**post_processing_params).execute()

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")

    input_file_type = 'fasta'
    base_output_dir = os.path.abspath(args.output_dir)

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = pvactools.lib.run_utils.split_algorithms(args.prediction_algorithms)
    (class_i_alleles, class_ii_alleles, species) = pvactools.lib.run_utils.split_alleles(args.allele)
    class_ii_alleles = pvactools.lib.run_utils.combine_class_ii_alleles(class_ii_alleles)

    shared_arguments = {
        'input_file'                : args.input_file,
        'input_file_type'           : input_file_type,
        'sample_name'               : args.sample_name,
        'top_score_metric'          : args.top_score_metric,
        'binding_threshold'         : args.binding_threshold,
        'percentile_threshold'      : args.percentile_threshold,
        'allele_specific_binding_thresholds': args.allele_specific_binding_thresholds,
        'net_chop_fasta'            : args.input_file,
        'net_chop_method'           : args.net_chop_method,
        'net_chop_threshold'        : args.net_chop_threshold,
        'additional_report_columns' : args.additional_report_columns,
        'fasta_size'                : args.fasta_size,
        'iedb_retries'              : args.iedb_retries,
        'keep_tmp_files'            : args.keep_tmp_files,
        'n_threads'                 : args.n_threads,
        'species'                   : species,
        'run_reference_proteome_similarity': args.run_reference_proteome_similarity,
        'blastp_path'               : args.blastp_path,
        'blastp_db'                 : args.blastp_db,
        'run_post_processor'        : True,
        'exclude_NAs'               : args.exclude_NAs,
    }

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0:
        if args.iedb_install_directory:
            iedb_mhc_i_executable = os.path.join(args.iedb_install_directory, 'mhc_i', 'src', 'predict_binding.py')
            if not os.path.exists(iedb_mhc_i_executable):
                sys.exit("IEDB MHC I executable path doesn't exist %s" % iedb_mhc_i_executable)
        else:
            iedb_mhc_i_executable = None

        print("Executing MHC Class I predictions")

        output_dir = os.path.join(base_output_dir, 'MHC_Class_I')
        os.makedirs(output_dir, exist_ok=True)

        class_i_arguments = shared_arguments.copy()
        class_i_arguments['alleles']                 = class_i_alleles
        class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
        class_i_arguments['epitope_lengths']         = args.class_i_epitope_length
        class_i_arguments['prediction_algorithms']   = class_i_prediction_algorithms
        class_i_arguments['output_dir']              = output_dir
        class_i_arguments['netmhc_stab']             = args.netmhc_stab
        pipeline = PvacbindPipeline(**class_i_arguments)
        pipeline.execute()
    elif len(class_i_prediction_algorithms) == 0:
        print("No MHC class I prediction algorithms chosen. Skipping MHC class I predictions.")
    elif len(class_i_alleles) == 0:
        print("No MHC class I alleles chosen. Skipping MHC class I predictions.")

    if len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        if args.iedb_install_directory:
            iedb_mhc_ii_executable = os.path.join(args.iedb_install_directory, 'mhc_ii', 'mhc_II_binding.py')
            if not os.path.exists(iedb_mhc_ii_executable):
                sys.exit("IEDB MHC II executable path doesn't exist %s" % iedb_mhc_ii_executable)
        else:
            iedb_mhc_ii_executable = None

        print("Executing MHC Class II predictions")

        output_dir = os.path.join(base_output_dir, 'MHC_Class_II')
        os.makedirs(output_dir, exist_ok=True)

        class_ii_arguments = shared_arguments.copy()
        class_ii_arguments['alleles']                 = class_ii_alleles
        class_ii_arguments['prediction_algorithms']   = class_ii_prediction_algorithms
        class_ii_arguments['iedb_executable']         = iedb_mhc_ii_executable
        class_ii_arguments['epitope_lengths']         = args.class_ii_epitope_length
        class_ii_arguments['output_dir']              = output_dir
        class_ii_arguments['netmhc_stab']             = False
        pipeline = PvacbindPipeline(**class_ii_arguments)
        pipeline.execute()
    elif len(class_ii_prediction_algorithms) == 0:
        print("No MHC class II prediction algorithms chosen. Skipping MHC class II predictions.")
    elif len(class_ii_alleles) == 0:
        print("No MHC class II alleles chosen. Skipping MHC class II predictions.")

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0 and len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        print("Creating combined reports")
        create_combined_reports(base_output_dir, args)

    pvactools.lib.run_utils.change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
