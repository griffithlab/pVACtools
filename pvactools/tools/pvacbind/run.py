import sys
import argparse
import os
import shutil
import yaml
import platform
import copy

from pvactools.lib.fasta_to_kmers import SequenceFastaToKmers
from pvactools.lib.prediction_class import *
from pvactools.lib.pipeline import PvacbindPipeline
from pvactools.lib.run_argument_parser import PvacbindRunArgumentParser
from pvactools.lib.post_processor import PostProcessor
from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *
from pvactools.lib.print_log import *

def define_parser():
    return PvacbindRunArgumentParser().parser

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
    post_processing_params['species'] = run_params['species']
    post_processing_params['file_type'] = 'pVACbind'
    PostProcessor(**post_processing_params).execute()

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")

    if args.n_threads > 1 and platform.system() == "Darwin":
        raise Exception("Multithreading is not supported on MacOS")

    if (args.netmhciipan_version == '4.0' and args.iedb_install_directory is not None):
        raise Exception("Standalone IEDB does not support version 4.0")
    NetMHCIIVersion.netmhciipan_version = args.netmhciipan_version

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    alleles = combine_class_ii_alleles(args.allele)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(alleles)

    input_file_type = 'fasta'
    base_output_dir = os.path.abspath(args.output_dir)
    os.makedirs(base_output_dir, exist_ok=True)

    print_log(os.path.join(base_output_dir, 'log'), vars(args), 'inputs')

    for length in set(args.class_i_epitope_length + args.class_ii_epitope_length):
        fasta_to_kmer_arguments = {
            'fasta': args.input_file,
            'output_dir': base_output_dir,
            'epitope_length': length,
            'sample_name': args.sample_name,
        }
        SequenceFastaToKmers(**fasta_to_kmer_arguments).execute()

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
        'problematic_amino_acids'   : args.problematic_amino_acids,
        'run_post_processor'        : True,
        'peptide_fasta'             : args.peptide_fasta,
        'aggregate_inclusion_binding_threshold': args.aggregate_inclusion_binding_threshold,
        'aggregate_inclusion_count_limit': args.aggregate_inclusion_count_limit,
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
                pipeline = PvacbindPipeline(**per_length_run_arguments)
                pipeline.execute()
                output_file = os.path.join(per_epitope_output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
                if os.path.exists(output_file):
                    output_files.append(output_file)
            if len(output_files) > 0:
                all_epitopes_file = os.path.join(output_dir, "{}.MHC_{}.all_epitopes.tsv".format(args.sample_name,mhc_class))
                filtered_file = os.path.join(output_dir, "{}.MHC_{}.filtered.tsv".format(args.sample_name,mhc_class))
                post_processing_params = vars(args).copy()
                post_processing_params['fasta'] = args.input_file
                post_processing_params["filename_addition"] = "MHC_{}".format(mhc_class)
                create_per_class_report(output_files, all_epitopes_file, filtered_file, post_processing_params, run_arguments)
            else:
                print("\nNo processable variants found. Aborting.\n")
        elif len(prediction_algorithms) == 0:
            print("No MHC class {} prediction algorithms chosen. Skipping MHC class {} predictions.".format(mhc_class, mhc_class))
        elif len(alleles) == 0:
            print("No MHC class {} alleles chosen. Skipping MHC class {} predictions.".format(mhc_class, mhc_class))

    change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
