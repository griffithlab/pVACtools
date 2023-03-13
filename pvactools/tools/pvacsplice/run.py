import sys
import os
import pandas as pd
from pathlib import Path

from pvactools.lib.splice_pipeline import *
from pvactools.lib.prediction_class import *
from pvactools.lib.pipeline import *
from pvactools.lib.run_argument_parser import *
from pvactools.lib.post_processor import *
from pvactools.lib.run_utils import *
from pvactools.lib.print_log import *


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
    combined_output_file = os.path.join(output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
    pvactools.lib.run_utils.combine_reports([file1, file2], combined_output_file)
    filtered_report_file = os.path.join(output_dir, "{}.filtered.tsv".format(args.sample_name))

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    # fasta
    if Path(args.ref_fasta).suffix not in ['.fa', '.fasta']:
        sys.exit('The fasta input path does not point to a fasta file.')
    elif Path(args.ref_fasta).suffix == '.gz':
        sys.exit('pVACsplice does not currently support gzipped fasta files.')
    # gtf
    if Path(args.gtf_file).suffix not in ['.gtf', '.gz']:
        sys.exit('The gtf input path does not point to a gtf file.')
    # vcf
    if Path(args.annotated_vcf).suffix not in ['.vcf', '.gz']:
        sys.exit('The vcf input path does not point to a vcf file.')
    # vcf gz.tbi index file
    if not Path(f'{args.annotated_vcf}.tbi'):
        sys.exit('Gzipped VCF files must be indexed. (tabix -p vcf <vcf_file>)')
    # iedb retries - default 5
    if args.iedb_retries > 100:
        sys.exit("The number of IEDB retries must be less than or equal to 100")
    # fasta size
    if args.fasta_size % 2 != 0:
        sys.exit("The fasta size needs to be an even number")

    base_output_dir = os.path.abspath(args.output_dir) # junctions dir
    os.makedirs(base_output_dir, exist_ok=True)

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(args.allele)

    # input file check
    print_log(os.path.join(base_output_dir, 'log'), vars(args))

    junction_arguments = {
        'input_file_type'                  : 'junctions',
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
        'anchor_types'                     : args.anchor_types,
        'save_gtf'                         : args.save_gtf,
        'normal_sample_name'               : args.normal_sample_name,
        'class_i_hla'                      : class_i_alleles,
        'class_ii_hla'                     : class_ii_alleles,
    }

    pipeline = JunctionPipeline(**junction_arguments)
    pipeline.execute()

    splice_arguments = junction_arguments.copy()
    additional_args = {
        'splice_output_dir'         : base_output_dir,
        'output_dir'                : os.path.join(base_output_dir, 'MHC_Class_I'), # todo remove 'I' and make compatible for both classes
        'top_score_metric'          : args.top_score_metric,
        'binding_threshold'         : args.binding_threshold,
        'percentile_threshold'      : args.percentile_threshold,
        'allele_specific_binding_thresholds': args.allele_specific_binding_thresholds,
        'aggregate_inclusion_binding_threshold' : args.aggregate_inclusion_binding_threshold,
        'net_chop_method'           : args.net_chop_method,
        'net_chop_threshold'        : args.net_chop_threshold,
        'additional_report_columns' : args.additional_report_columns,
        'fasta_size'                : args.fasta_size,
        'iedb_retries'              : args.iedb_retries,
        'keep_tmp_files'            : args.keep_tmp_files,
        'n_threads'                 : args.n_threads,
        'species'                   : species,
        'run_reference_proteome_similarity': args.run_reference_proteome_similarity,
        'blastp_db'                 : args.blast_dp, # todo add these to run_args_parser
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
        'run_post_processor'        : True,
        'exclude_NAs'               : args.exclude_NAs,
    }
    splice_arguments.update(additional_args)

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0:
        if args.iedb_install_directory:
            iedb_mhc_i_executable = os.path.join(args.iedb_install_directory, 'mhc_i', 'src', 'predict_binding.py')
            if not os.path.exists(iedb_mhc_i_executable):
                sys.exit("IEDB MHC I executable path doesn't exist %s" % iedb_mhc_i_executable)
        else:
            iedb_mhc_i_executable = None
        
        for x in args.class_i_epitope_length: # todo reformat loops here after comparing to pvacfuse run

            print(f'Executing MHC Class I predictions for {x}mers')
            output_dir = os.path.join(base_output_dir, 'MHC_Class_I', f'MHC_Class_I_{x}')
            os.makedirs(output_dir, exist_ok=True)

            class_i_arguments = splice_arguments.copy()
            class_i_arguments['input_file']              = f'{base_output_dir}/{args.sample_name}.{x}.fa'
            class_i_arguments['alleles']                 = class_i_alleles
            class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
            class_i_arguments['epitope_lengths']         = x
            class_i_arguments['prediction_algorithms']   = class_i_prediction_algorithms
            class_i_arguments['output_dir']              = output_dir
            class_i_arguments['netmhc_stab']             = args.netmhc_stab
            
            pipeline = PvacsplicePipeline(**class_i_arguments)
            pipeline.execute()

        combine_reports_per_class(base_output_dir, args, 'I', class_i_prediction_algorithms)

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

            class_ii_arguments = splice_arguments.copy()
            class_ii_arguments['input_file']              = f'{base_output_dir}/{args.sample_name}.{y}.fa'
            class_ii_arguments['alleles']                 = class_ii_alleles
            class_ii_arguments['prediction_algorithms']   = class_ii_prediction_algorithms
            class_ii_arguments['iedb_executable']         = iedb_mhc_ii_executable
            class_ii_arguments['epitope_lengths']         = y
            class_ii_arguments['output_dir']              = output_dir
            class_ii_arguments['netmhc_stab']             = False
            
            pipeline = PvacsplicePipeline(**class_ii_arguments)
            pipeline.execute()

        combine_reports_per_class(base_output_dir, args, 'II', class_ii_prediction_algorithms)

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
