import sys
import argparse
import os
import shutil
import yaml
import platform
import copy

from pvactools.lib.fasta_to_kmers import SequenceFastaToKmers
from pvactools.lib.pvacbind_prediction_pipeline import PvacbindPredictionPipeline
from pvactools.lib.prediction_class import *
from pvactools.lib.run_argument_parser import PvacbindRunArgumentParser
from pvactools.lib.post_processor import PostProcessor
from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *
from pvactools.lib.print_log import *

def define_parser():
    return PvacbindRunArgumentParser().parser

def create_per_class_report(files, all_epitopes_output_file, filtered_report_file, post_processing_params):
    for file_name in files:
        if not os.path.exists(file_name):
            print("File {} doesn't exist. Aborting.".format(file_name))
            return

    combine_reports(files, all_epitopes_output_file)

    post_processing_params['input_file'] = all_epitopes_output_file
    post_processing_params['filtered_report_file'] = filtered_report_file
    post_processing_params['run_coverage_filter'] = False
    post_processing_params['run_transcript_support_level_filter'] = False
    post_processing_params['run_manufacturability_metrics'] = True
    post_processing_params['run_net_chop'] = True if post_processing_params['net_chop_method'] else False
    post_processing_params['run_netmhc_stab'] = True if post_processing_params['netmhc_stab'] else False
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
        if len(params['prediction_algorithms']) > 0 and len(params['alleles']) > 0:
            params['base_output_dir'] = base_output_dir
            params['mhc_class'] = mhc_class
            params['sample_name'] = args.sample_name
            params['fasta_size'] = args.fasta_size
            params['iedb_retries'] = args.iedb_retries
            params['n_threads'] = args.n_threads
            params['additional_report_columns'] = args.additional_report_columns

            predictor = PvacbindPredictionPipeline(**params)
            predictor.execute()

            if len(predictor.output_files) > 0:
                all_epitopes_file = os.path.join(predictor.output_dir, "{}.MHC_{}.all_epitopes.tsv".format(args.sample_name,mhc_class))
                filtered_file = os.path.join(predictor.output_dir, "{}.MHC_{}.filtered.tsv".format(args.sample_name,mhc_class))
                post_processing_params = vars(args).copy()
                post_processing_params['fasta'] = args.input_file
                post_processing_params['net_chop_fasta'] = args.input_file
                post_processing_params["filename_addition"] = "MHC_{}".format(mhc_class)
                post_processing_params["species"] = species
                post_processing_params["netmhc_stab"] = params["netmhc_stab"]
                create_per_class_report(predictor.output_files, all_epitopes_file, filtered_file, post_processing_params)
            else:
                print("\nNo processable variants found. Aborting.\n")
        elif len(params["prediction_algorithms"]) == 0:
            print("No MHC class {} prediction algorithms chosen. Skipping MHC class {} predictions.".format(mhc_class, mhc_class))
        elif len(params["alleles"]) == 0:
            print("No MHC class {} alleles chosen. Skipping MHC class {} predictions.".format(mhc_class, mhc_class))

    change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
