import sys
import argparse
import os
import shutil
import platform
import copy

from pvactools.lib.prediction_class import *
from pvactools.lib.pvacfuse_prediction_pipeline import PvacfusePredictionPipeline
from pvactools.lib.fusion_to_kmer_pipeline import FusionToKmerPipeline
from pvactools.lib.run_argument_parser import PvacfuseRunArgumentParser
from pvactools.lib.post_processor import PostProcessor
from pvactools.tools.pvacfuse.generate_protein_fasta import PvacfuseGenerateProteinFasta
from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *
from pvactools.lib.print_log import *

def define_parser():
    return PvacfuseRunArgumentParser().parser

def create_per_class_report(files, all_epitopes_output_file, filtered_report_file, post_processing_params):
    for file_name in files:
        if not os.path.exists(file_name):
            print("File {} doesn't exist. Aborting.".format(file_name))
            return

    combine_reports(files, all_epitopes_output_file)

    post_processing_params['input_file'] = all_epitopes_output_file
    post_processing_params['filtered_report_file'] = filtered_report_file
    post_processing_params['run_coverage_filter'] = True
    post_processing_params['run_transcript_support_level_filter'] = False
    post_processing_params['run_manufacturability_metrics'] = True
    post_processing_params['run_net_chop'] = True if post_processing_params['net_chop_method'] else False
    post_processing_params['run_netmhc_stab'] = True if post_processing_params['netmhc_stab'] else False
    post_processing_params['file_type'] = 'pVACfuse'
    PostProcessor(**post_processing_params).execute()

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

    if args.n_threads > 1 and platform.system() == "Darwin":
        raise Exception("Multithreading is not supported on MacOS")

    base_output_dir = os.path.abspath(args.output_dir)

    if (args.netmhciipan_version == '4.0' and args.iedb_install_directory is not None):
        raise Exception("Standalone IEDB does not support version 4.0")
    NetMHCIIVersion.netmhciipan_version = args.netmhciipan_version

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    alleles = combine_class_ii_alleles(args.allele)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(alleles)

    fusions_dir = os.path.abspath(args.output_dir)
    os.makedirs(fusions_dir, exist_ok=True)

    print_log(os.path.join(fusions_dir, 'log'), vars(args), 'inputs')

    fusion_arguments = {
        'output_dir'              : fusions_dir,
        'input_file'              : args.input_file,
        'sample_name'             : args.sample_name,
        'transcript_fasta'        : args.ref_fasta,
        'starfusion_file'         : args.starfusion_file,
        'class_i_epitope_length'  : args.class_i_epitope_length,
        'class_ii_epitope_length' : args.class_ii_epitope_length,
        'class_i_hla'             : class_i_alleles,
        'class_ii_hla'            : class_ii_alleles,
    }

    pipeline = FusionToKmerPipeline(**fusion_arguments)
    pipeline.execute()

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
            params['transcript_fasta'] = pipeline.create_file_path('fasta')

            predictor = PvacfusePredictionPipeline(**params)
            predictor.execute()

            if len(predictor.output_files) > 0:
                post_processing_params = vars(args).copy()
                if args.run_reference_proteome_similarity:
                    fasta_file = os.path.join(predictor.output_dir, "{}.7.fasta".format(args.sample_name))
                    if not os.path.exists(fasta_file):
                        fasta_params = {
                            'fasta_file_path': pipeline.create_file_path('fasta'),
                            'trimmed_fasta_file_path': fasta_file,
                            'flanking_sequence_length': 7,
                            'mutant_only': False,
                        }
                        PvacfuseGenerateProteinFasta(**fasta_params).trim_sequences()
                    post_processing_params['fasta'] = fasta_file
                # generate and copy net_chop fasta to output dir if specified
                if args.net_chop_method:
                    fasta_file = os.path.join(predictor.output_dir, "{}.10.fasta".format(args.sample_name))
                    if not os.path.exists(fasta_file):
                        fasta_params = {
                            'fasta_file_path': pipeline.create_file_path('fasta'),
                            'trimmed_fasta_file_path': fasta_file,
                            'flanking_sequence_length': 10,
                            'mutant_only': False,
                        }
                        PvacfuseGenerateProteinFasta(**fasta_params).trim_sequences()
                    post_processing_params['net_chop_fasta'] = fasta_file
                all_epitopes_file = os.path.join(predictor.output_dir, "{}.MHC_{}.all_epitopes.tsv".format(args.sample_name,mhc_class))
                filtered_file = os.path.join(predictor.output_dir, "{}.MHC_{}.filtered.tsv".format(args.sample_name,mhc_class))
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
