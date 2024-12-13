import sys
import os
import pandas as pd
import shutil
from pathlib import Path
from pvactools.lib.splice_pipeline import *
from pvactools.lib.prediction_class import *
from pvactools.lib.pipeline import *
from pvactools.lib.run_argument_parser import *
from pvactools.lib.post_processor import *
from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *
from pvactools.lib.print_log import *


def define_parser():
    return PvacspliceRunArgumentParser().parser

def combine_epitope_len_reports(file_list, file_final_name):
    if len(file_list) > 1:
        combine_reports(file_list, file_final_name)
    elif len(file_list) == 1:
        shutil.copy(file_list[0], file_final_name)

def create_full_combined_reports(base_output_dir, args):
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
    post_processing_params['minimum_fold_change'] = None
    post_processing_params['run_coverage_filter'] = True
    post_processing_params['run_transcript_support_level_filter'] = True
    post_processing_params['run_net_chop'] = False
    post_processing_params['run_netmhc_stab'] = False
    post_processing_params['run_manufacturability_metrics'] = False
    post_processing_params['run_reference_proteome_similarity'] = False
    post_processing_params['file_type'] = 'pVACsplice'

    PostProcessor(**post_processing_params).execute()

def combine_reports_per_class(class_output_dir:str, params:dict, mhc_class:str):
    output_dir = os.path.join(class_output_dir, f'MHC_Class_{mhc_class}')

    for x in ['all_epitopes', 'filtered']:
        mhc_dirs = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.startswith(f'MHC_Class_{mhc_class}')]
    if not mhc_dirs:
        print(f'MHC_Class_{mhc_class} subfolder(s) are missing')
    combined_files = [os.path.join(m, f'{params["sample_name"]}.all_epitopes.tsv') for m in mhc_dirs]
    combined_fn = os.path.join(output_dir, f'{params["sample_name"]}.all_epitopes.tsv')
    combine_epitope_len_reports(combined_files, combined_fn)
    filtered_fn = os.path.join(output_dir, f'{params["sample_name"]}.filtered.tsv')

    post_processing_params = params.copy()
    post_processing_params['file_type'] = 'pVACsplice'
    post_processing_params['input_file'] = combined_fn
    post_processing_params['filtered_report_file'] = filtered_fn
    post_processing_params['minimum_fold_change'] = None
    # methods in pp class
    post_processing_params['run_manufacturability_metrics'] = True
    post_processing_params['run_coverage_filter'] = True
    post_processing_params['run_transcript_support_level_filter'] = True
    # add custom params for netchop / netmhc_stab
    if params['net_chop_method']:
        post_processing_params['net_chop_fasta'] = params['net_chop_fasta']
        post_processing_params['run_net_chop'] = True
    else:
        post_processing_params['run_net_chop'] = False
    if params['netmhc_stab']:
        post_processing_params['run_netmhc_stab'] = True
    else:
        post_processing_params['run_netmhc_stab'] = False

    print('Begin post processor')
    PostProcessor(**post_processing_params).execute()

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    # fasta
    if Path(args.ref_fasta).suffix not in ['.fa', '.fasta']:
        sys.exit('The fasta input path does not point to a fasta file.')
    if is_gz_file(Path(args.ref_fasta)):
        sys.exit('pVACsplice does not currently support gzipped fasta files.')
    # supplied fasta size length should be even
    if args.fasta_size%2 != 0:
        sys.exit("The fasta size needs to be an even number")
    # gtf
    if Path(args.gtf_file).suffix not in ['.gtf', '.tsv'] and not is_gz_file(args.gtf_file):
        sys.exit('The gtf input path does not point to a gtf file.')
    # vcf
    if Path(args.annotated_vcf).suffix != '.vcf' and not is_gz_file(args.annotated_vcf):
        sys.exit('The vcf input path does not point to a vcf file.')
    # vcf gz.tbi index file
    if is_gz_file(args.annotated_vcf) and not Path(f'{args.annotated_vcf}.tbi'):
        sys.exit('Gzipped VCF files must be indexed. (tabix -p vcf <vcf_file>)')
    # iedb retries - default 5
    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100.")

    # pvacsplice output dir (from args)
    junctions_dir = os.path.abspath(args.output_dir)
    os.makedirs(junctions_dir, exist_ok=True)

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(args.allele)

    # all input file check
    print_log(os.path.join(junctions_dir, 'log'), vars(args), 'inputs')

    junction_arguments = {
        'input_file_type'                  : 'junctions',
        'junctions_dir'                    : junctions_dir,
        'input_file'                       : args.input_file,
        'gtf_file'                         : args.gtf_file,
        'save_gtf'                         : args.save_gtf,
        'sample_name'                      : args.sample_name,
        'ref_fasta'                        : args.ref_fasta,
        'annotated_vcf'                    : args.annotated_vcf,
        'pass_only'                        : args.pass_only,
        'class_i_epitope_length'           : args.class_i_epitope_length,
        'class_ii_epitope_length'          : args.class_ii_epitope_length,
        'biotypes'                         : args.biotypes,
        'junction_score'                   : args.junction_score,
        'variant_distance'                 : args.variant_distance,
        'anchor_types'                     : args.anchor_types,
        'normal_sample_name'               : args.normal_sample_name,
        'class_i_hla'                      : class_i_alleles,
        'class_ii_hla'                     : class_ii_alleles,
        'keep_tmp_files'                   : args.keep_tmp_files,
    }

    pipeline = JunctionPipeline(**junction_arguments)
    pipeline.execute()


    additional_args = {
        'top_score_metric'          : args.top_score_metric,
        'binding_threshold'         : args.binding_threshold,
        'percentile_threshold'      : args.percentile_threshold,
        'allele_specific_binding_thresholds': args.allele_specific_binding_thresholds,
        'aggregate_inclusion_binding_threshold' : args.aggregate_inclusion_binding_threshold,
        'aggregate_inclusion_count_limit': args.aggregate_inclusion_count_limit,
        'net_chop_method'           : args.net_chop_method,
        'net_chop_threshold'        : args.net_chop_threshold,
        'additional_report_columns' : args.additional_report_columns,
        'fasta_size'                : args.fasta_size,
        'iedb_retries'              : args.iedb_retries,
        'n_threads'                 : args.n_threads,
        'species'                   : species,
        'run_reference_proteome_similarity': args.run_reference_proteome_similarity,
        'blastp_db'                 : args.blastp_db,
        'blastp_path'               : args.blastp_path,
        'peptide_fasta'             : args.peptide_fasta,
        'problematic_amino_acids'   : args.problematic_amino_acids,
        'normal_cov'                : args.normal_cov,
        'normal_vaf'                : args.normal_vaf,
        'tdna_cov'                  : args.tdna_cov,
        'tdna_vaf'                  : args.tdna_vaf,
        'trna_cov'                  : args.trna_cov,
        'trna_vaf'                  : args.trna_vaf,
        'expn_val'                  : args.expn_val,
        'tumor_purity'              : args.tumor_purity,
        'maximum_transcript_support_level' : args.maximum_transcript_support_level,
        'run_post_processor'        : True,
        'exclude_NAs'               : args.exclude_NAs,
    }
    junction_arguments.update(additional_args)

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0:
        if args.iedb_install_directory:
            iedb_mhc_i_executable = os.path.join(args.iedb_install_directory, 'mhc_i', 'src', 'predict_binding.py')
            if not os.path.exists(iedb_mhc_i_executable):
                sys.exit("IEDB MHC I executable path doesn't exist %s" % iedb_mhc_i_executable)
        else:
            iedb_mhc_i_executable = None

        for x in args.class_i_epitope_length:

            print(f'Executing MHC Class I predictions for {x}mers')
            output_len_dir = os.path.join(junctions_dir, 'MHC_Class_I', f'MHC_Class_I_{x}')
            os.makedirs(output_len_dir, exist_ok=True)

            class_i_arguments = junction_arguments.copy()
            class_i_arguments['input_file']              = os.path.join(junctions_dir, 'tmp', f'{args.sample_name}.{x}.fa')
            class_i_arguments['alleles']                 = class_i_alleles
            class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
            class_i_arguments['epitope_lengths']         = x
            class_i_arguments['prediction_algorithms']   = class_i_prediction_algorithms
            class_i_arguments['output_dir']              = output_len_dir
            class_i_arguments['netmhc_stab']             = args.netmhc_stab

            pipeline = PvacsplicePipeline(**class_i_arguments)
            pipeline.execute()

        fasta_file = os.path.join(junctions_dir, "{}.transcripts.fa".format(args.sample_name))
        class_i_arguments['fasta'] = fasta_file
        class_i_arguments['net_chop_fasta'] = fasta_file

        combine_reports_per_class(junctions_dir, class_i_arguments, 'I')

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
            output_len_dir = os.path.join(junctions_dir, 'MHC_Class_II', f'MHC_Class_II_{y}')
            os.makedirs(output_len_dir, exist_ok=True)

            class_ii_arguments = junction_arguments.copy()
            class_ii_arguments['input_file']              = os.path.join(junctions_dir, 'tmp', f'{args.sample_name}.{y}.fa')
            class_ii_arguments['alleles']                 = class_ii_alleles
            class_ii_arguments['prediction_algorithms']   = class_ii_prediction_algorithms
            class_ii_arguments['iedb_executable']         = iedb_mhc_ii_executable
            class_ii_arguments['epitope_lengths']         = y
            class_ii_arguments['output_dir']              = output_len_dir
            class_ii_arguments['netmhc_stab']             = False

            pipeline = PvacsplicePipeline(**class_ii_arguments)
            pipeline.execute()

        fasta_file = os.path.join(junctions_dir, "{}.transcripts.fa".format(args.sample_name))
        class_ii_arguments['fasta'] = fasta_file
        class_ii_arguments['net_chop_fasta'] = fasta_file

        combine_reports_per_class(junctions_dir, class_ii_arguments, 'II')

    elif len(class_ii_prediction_algorithms) == 0:
        print("No MHC class II prediction algorithms chosen. Skipping MHC class II predictions.")
    elif len(class_ii_alleles) == 0:
        print("No MHC class II alleles chosen. Skipping MHC class II predictions.")

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0 and len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        print("Creating combined reports")
        create_full_combined_reports(junctions_dir, args)

    if args.save_gtf is False:
        shutil.rmtree(os.path.join(junctions_dir, 'tmp'), ignore_errors=True)

    change_permissions_recursive(junctions_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
