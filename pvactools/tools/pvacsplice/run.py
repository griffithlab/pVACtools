import sys
import os
import pandas as pd
from pvactools.tools.pvacsplice.splice_pipeline import *
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

    file1 = os.path.join(base_output_dir, 'MHC_Class_I', "{}.all_epitopes.tsv".format(args.sample_name))
    file2 = os.path.join(base_output_dir, 'MHC_Class_II', "{}.all_epitopes.tsv".format(args.sample_name))
    if not os.path.exists(file1):
        print("File {} doesn't exist. Aborting.".format(file1))
        return
    if not os.path.exists(file2):
        print("File {} doesn't exist. Aborting.".format(file2))
        return
    combined_output_file = os.path.join(output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
    combine_reports([file1, file2], combined_output_file)
    filtered_report_file = os.path.join(output_dir, "{}.filtered.tsv".format(args.sample_name))

    post_processing_params = vars(args)
    post_processing_params['input_file'] = combined_output_file
    post_processing_params['filtered_report_file'] = filtered_report_file
    post_processing_params['run_coverage_filter'] = False
    post_processing_params['minimum_fold_change'] = None
    post_processing_params['file_type'] = 'pVACsplice'
    post_processing_params['run_transcript_support_level_filter'] = False
    post_processing_params['run_net_chop'] = False
    post_processing_params['run_netmhc_stab'] = False
    post_processing_params['run_manufacturability_metrics'] = False
    post_processing_params['run_reference_proteome_similarity'] = False

    PostProcessor(**post_processing_params).execute()

def combine_reports_mhc_class(base_output_dir, args, mhc_class):
    output_dir = os.path.join(base_output_dir, mhc_class)
    
    files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith('all_epitopes.tsv')]
    print(files)
    combined_file = os.path.join(output_dir, f'{args.sample_name}.all_epitopes.tsv')
    filtered_file = os.path.join(output_dir, f'{args.sample_name}.filtered.tsv')
    
    if len(files) > 1:
        combine_reports(files, combined_file)
        print('Finish combined report')
        if os.path.exists(combined_file):
            for f in files:
                os.remove(f)
    elif len(files) == 1:
        os.rename(files[0], combined_file)

    post_processing_params = vars(args)
    post_processing_params['input_file'] = combined_file
    post_processing_params['file_type'] = 'pVACsplice'
    post_processing_params['filtered_report_file'] = filtered_file
    post_processing_params['run_coverage_filter'] = True
    post_processing_params['run_transcript_support_level_filter'] = False
    post_processing_params['minimum_fold_change'] = None
    post_processing_params['run_manufacturability_metrics'] = True
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

    base_output_dir = os.path.abspath(args.output_dir)
    os.makedirs(base_output_dir, exist_ok=True)
    tmp_dir = os.path.join(base_output_dir, 'tmp')

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(args.allele)


    junction_arguments = {
        'input_file'                       : args.input_file, 
        'sample_name'                      : args.sample_name,
        'base_output_dir'                  : base_output_dir,
        'ref_fasta'                        : args.ref_fasta,
        'annotated_vcf'                    : args.annotated_vcf,
        'class_i_epitope_length'           : args.class_i_epitope_length,
        'class_ii_epitope_length'          : args.class_ii_epitope_length,
        'maximum_transcript_support_level' : args.maximum_transcript_support_level,
        'junction_score'                   : args.junction_score,
        'variant_distance'                 : args.variant_distance,
        'normal_sample_name'               : args.normal_sample_name,
    }

    pipeline = JunctionPipeline(**junction_arguments)
    pipeline.execute()

    pvacbind_arguments = junction_arguments.copy()
    additional_args = {
        'input_file_type'           : 'junctions',
        'base_output_dir'           : base_output_dir,
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
    pvacbind_arguments.update(additional_args)

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

        for x in args.class_i_epitope_length:

            class_i_arguments = pvacbind_arguments.copy()
            class_i_arguments['input_file']              = f'{tmp_dir}/peptides_length_{x}.fa'
            class_i_arguments['alleles']                 = class_i_alleles
            class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
            class_i_arguments['epitope_lengths']         = x
            class_i_arguments['prediction_algorithms']   = class_i_prediction_algorithms
            class_i_arguments['output_dir']              = output_dir
            class_i_arguments['netmhc_stab']             = args.netmhc_stab
            
            pipeline = PvacsplicePipeline(**class_i_arguments)
            pipeline.execute()

        combine_reports_mhc_class(base_output_dir, args, 'MHC_Class_I')

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

        for y in args.class_ii_epitope_length:

            class_ii_arguments = pvacbind_arguments.copy()
            class_i_arguments['input_file']               = f'{tmp_dir}/peptides_length_{y}.fa'
            class_ii_arguments['alleles']                 = class_ii_alleles
            class_ii_arguments['prediction_algorithms']   = class_ii_prediction_algorithms
            class_ii_arguments['iedb_executable']         = iedb_mhc_ii_executable
            class_ii_arguments['epitope_lengths']         = y
            class_ii_arguments['output_dir']              = output_dir
            class_ii_arguments['netmhc_stab']             = False
            
            pipeline = PvacsplicePipeline(**class_ii_arguments)
            pipeline.execute()

        combine_reports_mhc_class(base_output_dir, args, 'MHC_Class_II')

    elif len(class_ii_prediction_algorithms) == 0:
        print("No MHC class II prediction algorithms chosen. Skipping MHC class II predictions.")
    elif len(class_ii_alleles) == 0:
        print("No MHC class II alleles chosen. Skipping MHC class II predictions.")


    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0 and len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        print("Creating combined reports")
        create_combined_reports(base_output_dir, args)

    change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
