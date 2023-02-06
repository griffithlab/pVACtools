import sys
import os
import pandas as pd
from pvactools.lib.splice_pipeline import *
from pvactools.lib.prediction_class import *
from pvactools.lib.pipeline import *
from pvactools.lib.run_argument_parser import *
from pvactools.lib.post_processor import *
from pvactools.lib.run_utils import *


def define_parser():
    return PvacspliceRunArgumentParser().parser

def create_combined_reports(base_output_dir, args):
    output_dir = os.path.join(base_output_dir, 'combined')
    os.makedirs(output_dir, exist_ok=True)

    for x in ['all_epitopes', 'filtered']:
        file1 = os.path.join(base_output_dir, 'MHC_Class_I', f"{args.sample_name}.{x}.tsv")
        file2 = os.path.join(base_output_dir, 'MHC_Class_II', f"{args.sample_name}.{x}.tsv")
        if not os.path.exists(file1):
            print("File {} doesn't exist. Aborting.".format(file1))
            return
        if not os.path.exists(file2):
            print("File {} doesn't exist. Aborting.".format(file2))
            return
        output_file = os.path.join(output_dir, f"{args.sample_name}.{x}.tsv")
        combine_reports([file1, file2], output_file)


def combine_class_reports(file_list, file_final_name):
    if len(file_list) > 1:
        combine_reports(file_list, file_final_name)
    elif len(file_list) == 1:
        os.rename(file_list[0], file_final_name)


def combine_reports_per_class(base_output_dir, args, mhc_class):
    output_dir = os.path.join(base_output_dir, f'MHC_Class_{mhc_class}')

    mhc_dirs = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.startswith(f'MHC_Class_{mhc_class}')]

    if not mhc_dirs:
        print(f'MHC_Class_{mhc_class} subfolder(s) are missing')
    combined_files = [os.path.join(m, f'{args.sample_name}.all_epitopes.tsv') for m in mhc_dirs]

    combined_name = os.path.join(output_dir, f'{args.sample_name}.all_epitopes.tsv')
    filtered_name = os.path.join(output_dir, f'{args.sample_name}.filtered.tsv')

    combine_class_reports(combined_files, combined_name)

    post_processing_params = vars(args)
    post_processing_params['input_file'] = combined_name
    post_processing_params['file_type'] = 'pVACsplice'
    post_processing_params['filtered_report_file'] = filtered_name
    post_processing_params['run_coverage_filter'] = True
    post_processing_params['run_transcript_support_level_filter'] = False
    post_processing_params['minimum_fold_change'] = None
    post_processing_params['run_manufacturability_metrics'] = True
    post_processing_params['run_reference_proteome_similarity'] = False
    if args.net_chop_method:
        post_processing_params['net_chop_fasta'] = args.net_chop_fasta
        post_processing_params['run_net_chop'] = True
    else:
        post_processing_params['run_net_chop'] = False
    if args.netmhc_stab:
        post_processing_params['run_netmhc_stab'] = True
    else:
        post_processing_params['run_netmhc_stab'] = False
        
    print('Begin post processor')
    PostProcessor(**post_processing_params).execute()

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")

    base_output_dir = os.path.abspath(args.output_dir) # junctions dir
    os.makedirs(base_output_dir, exist_ok=True)

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(args.allele)

    junction_arguments = {
        'input_file'                       : args.input_file, 
        'gtf_file'                         : args.gtf_file,
        'sample_name'                      : args.sample_name,
        'base_output_dir'                  : base_output_dir,
        'ref_fasta'                        : args.ref_fasta,
        'annotated_vcf'                    : args.annotated_vcf,
        'class_i_epitope_length'           : args.class_i_epitope_length,
        'class_ii_epitope_length'          : args.class_ii_epitope_length,
        'maximum_transcript_support_level' : args.maximum_transcript_support_level,
        'junction_score'                   : args.junction_score,
        'variant_distance'                 : args.variant_distance,
        'save_gtf'                         : args.save_gtf,
        'normal_sample_name'               : args.normal_sample_name,
        'class_i_hla'                      : class_i_alleles,
        'class_ii_hla'                     : class_ii_alleles,
    }

    pipeline = JunctionPipeline(**junction_arguments)
    pipeline.execute()

    pvacsplice_arguments = junction_arguments.copy()
    additional_args = {
        'input_file_type'           : 'junctions',
        'splice_output_dir'         : base_output_dir,
        'output_dir'                : os.path.join(base_output_dir, 'MHC_Class_I'),
        'top_score_metric'          : args.top_score_metric,
        'binding_threshold'         : args.binding_threshold,
        'percentile_threshold'      : args.percentile_threshold,
        'allele_specific_cutoffs'   : args.allele_specific_binding_thresholds,
        'net_chop_method'           : args.net_chop_method,
        'net_chop_threshold'        : args.net_chop_threshold,
        'additional_report_columns' : args.additional_report_columns,
        'fasta_size'                : args.fasta_size,
        'iedb_retries'              : args.iedb_retries,
        'keep_tmp_files'            : args.keep_tmp_files,
        'n_threads'                 : args.n_threads,
        'species'                   : species,
        'run_reference_proteome_similarity': args.run_reference_proteome_similarity,
        'normal_cov'                : args.normal_cov,
        'normal_vaf'                : args.normal_vaf,
        'tdna_cov'                  : args.tdna_cov,
        'tdna_vaf'                  : args.tdna_vaf,
        'trna_cov'                  : args.trna_cov,
        'trna_vaf'                  : args.trna_vaf,
        'expn_val'                  : args.expn_val,
    }
    pvacsplice_arguments.update(additional_args)

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0:
        if args.iedb_install_directory:
            iedb_mhc_i_executable = os.path.join(args.iedb_install_directory, 'mhc_i', 'src', 'predict_binding.py')
            if not os.path.exists(iedb_mhc_i_executable):
                sys.exit("IEDB MHC I executable path doesn't exist %s" % iedb_mhc_i_executable)
        else:
            iedb_mhc_i_executable = None
        
        for x in args.class_i_epitope_length:

            print(f'Executing MHC Class I predictions for {x}mers')
            output_dir = os.path.join(base_output_dir, 'MHC_Class_I', f'MHC_Class_I_{x}')
            os.makedirs(output_dir, exist_ok=True)

            class_i_arguments = pvacsplice_arguments.copy()
            class_i_arguments['input_file']              = f'{base_output_dir}/{args.sample_name}.{x}.fa'
            class_i_arguments['alleles']                 = class_i_alleles
            class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
            class_i_arguments['epitope_lengths']         = x
            class_i_arguments['prediction_algorithms']   = class_i_prediction_algorithms
            class_i_arguments['output_dir']              = output_dir
            class_i_arguments['netmhc_stab']             = args.netmhc_stab
            
            pipeline = PvacsplicePipeline(**class_i_arguments)
            pipeline.execute()

        combine_reports_per_class(base_output_dir, args, 'I')

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

        for y in args.class_ii_epitope_length:

            print(f'Executing MHC Class II predictions for {y}mers')
            output_dir = os.path.join(base_output_dir, 'MHC_Class_II', f'MHC_Class_II_{y}')
            os.makedirs(output_dir, exist_ok=True)

            class_ii_arguments = pvacsplice_arguments.copy()
            class_ii_arguments['input_file']              = f'{base_output_dir}/{args.sample_name}.{y}.fa'
            class_ii_arguments['alleles']                 = class_ii_alleles
            class_ii_arguments['prediction_algorithms']   = class_ii_prediction_algorithms
            class_ii_arguments['iedb_executable']         = iedb_mhc_ii_executable
            class_ii_arguments['epitope_lengths']         = y
            class_ii_arguments['output_dir']              = output_dir
            class_ii_arguments['netmhc_stab']             = False
            
            pipeline = PvacsplicePipeline(**class_ii_arguments)
            pipeline.execute()

        combine_reports_per_class(base_output_dir, args, 'II')

    elif len(class_ii_prediction_algorithms) == 0:
        print("No MHC class II prediction algorithms chosen. Skipping MHC class II predictions.")
    elif len(class_ii_alleles) == 0:
        print("No MHC class II alleles chosen. Skipping MHC class II predictions.")


    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0 and len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        print("Creating combined reports")
        combine_class_reports(base_output_dir, args)

    change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
