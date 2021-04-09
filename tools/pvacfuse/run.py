import sys
import argparse
import os
import shutil
from lib.prediction_class import *
from lib.pipeline import *
from lib.run_argument_parser import *
from lib.post_processor import *
from lib.run_utils import *
import tools.pvacfuse.generate_protein_fasta
import lib.call_iedb

def define_parser():
    return PvacfuseRunArgumentParser().parser

def create_combined_reports(files, all_epitopes_output_file, filtered_report_file, run_manufacturability_metrics, args):
    for file_name in files:
        if not os.path.exists(file_name):
            print("File {} doesn't exist. Aborting.".format(file_name))
            return

    combine_reports(files, all_epitopes_output_file)

    post_processing_params = vars(args).copy()
    post_processing_params['input_file'] = all_epitopes_output_file
    post_processing_params['filtered_report_file'] = filtered_report_file
    post_processing_params['minimum_fold_change'] = None
    post_processing_params['run_coverage_filter'] = False
    post_processing_params['run_transcript_support_level_filter'] = False
    post_processing_params['run_net_chop'] = False
    post_processing_params['run_netmhc_stab'] = False
    post_processing_params['run_manufacturability_metrics'] = run_manufacturability_metrics
    post_processing_params['run_reference_proteome_similarity'] = False
    post_processing_params['file_type'] = 'pVACfuse'

    PostProcessor(**post_processing_params).execute()

def generate_fasta(args, output_dir, epitope_length):
    per_epitope_output_dir = os.path.join(output_dir, str(epitope_length))
    os.makedirs(per_epitope_output_dir, exist_ok=True)
    output_file = os.path.join(per_epitope_output_dir, "{}.fa".format(args.sample_name))
    params = [
        args.input_file,
        str(epitope_length - 1),
        output_file,
    ]
    if args.downstream_sequence_length is not None:
        params.extend(["-d", str(args.downstream_sequence_length)])
    else:
        params.extend(["-d", 'full'])
    tools.pvacfuse.generate_protein_fasta.main(params, save_tsv_file=True)
    os.unlink("{}.manufacturability.tsv".format(output_file))
    return (output_file, per_epitope_output_dir)

def append_columns(intermediate_output_file, tsv_file, output_file):
    tsv_entries = {}
    with open(tsv_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh, delimiter="\t")
        for line in reader:
            tsv_entries[line['index']] = line

    with open(intermediate_output_file, 'r') as input_fh, open(output_file, 'w') as output_fh:
        reader = csv.DictReader(input_fh, delimiter="\t")
        fieldnames = ['Chromosome', 'Start', 'Stop', 'Transcript', 'Gene Name', 'Variant Type'] + reader.fieldnames
        writer = csv.DictWriter(output_fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for line in reader:
            matching_line = tsv_entries[line['Mutation']]
            line['Chromosome'] = matching_line['chromosome_name']
            line['Start'] = matching_line['start']
            line['Stop'] = matching_line['stop']
            line['Transcript'] = matching_line['transcript_name']
            line['Gene Name'] = matching_line['gene_name']
            line['Variant Type'] = matching_line['variant_type']
            writer.writerow(line)

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


    base_output_dir = os.path.abspath(args.output_dir)

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = split_algorithms(args.prediction_algorithms)
    (class_i_alleles, class_ii_alleles, species) = split_alleles(args.allele)

    shared_arguments = {
        'input_file_type'           : 'fasta',
        'sample_name'               : args.sample_name,
        'top_score_metric'          : args.top_score_metric,
        'binding_threshold'         : args.binding_threshold,
        'percentile_threshold'      : args.percentile_threshold,
        'allele_specific_cutoffs'   : args.allele_specific_binding_thresholds,
        'net_chop_method'           : args.net_chop_method,
        'net_chop_threshold'        : args.net_chop_threshold,
        'additional_report_columns' : args.additional_report_columns,
        'fasta_size'                : args.fasta_size,
        'iedb_retries'              : args.iedb_retries,
        'downstream_sequence_length': downstream_sequence_length,
        'keep_tmp_files'            : args.keep_tmp_files,
        'n_threads'                 : args.n_threads,
        'species'                   : species,
        'run_reference_proteome_similarity': args.run_reference_proteome_similarity,
        'run_post_processor'        : False
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

        output_files = []
        for epitope_length in args.class_i_epitope_length:
            (input_file, per_epitope_output_dir) = generate_fasta(args, output_dir, epitope_length)

            class_i_arguments = shared_arguments.copy()
            class_i_arguments['input_file']              = input_file
            class_i_arguments['alleles']                 = class_i_alleles
            class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
            class_i_arguments['epitope_lengths']         = [epitope_length]
            class_i_arguments['prediction_algorithms']   = class_i_prediction_algorithms
            class_i_arguments['output_dir']              = per_epitope_output_dir
            class_i_arguments['netmhc_stab']             = args.netmhc_stab
            pipeline = PvacbindPipeline(**class_i_arguments)
            pipeline.execute()
            intermediate_output_file = os.path.join(per_epitope_output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
            output_file = os.path.join(per_epitope_output_dir, "{}.all_epitopes.final.tsv".format(args.sample_name))
            append_columns(intermediate_output_file, "{}.tsv".format(input_file), output_file)
            output_files.append(output_file)
            if epitope_length == max(args.class_i_epitope_length):
                fasta_file = os.path.join(output_dir, "{}.fasta".format(args.sample_name))
                shutil.copy(input_file, fasta_file)
        all_epitopes_file = os.path.join(output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
        filtered_file = os.path.join(output_dir, "{}.filtered.tsv".format(args.sample_name))
        create_combined_reports(output_files, all_epitopes_file, filtered_file, True, args)
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

        output_files = []
        for epitope_length in args.class_ii_epitope_length:
            (input_file, per_epitope_output_dir) = generate_fasta(args, output_dir, epitope_length)

            class_ii_arguments = shared_arguments.copy()
            class_ii_arguments['input_file']              = input_file
            class_ii_arguments['alleles']                 = class_ii_alleles
            class_ii_arguments['iedb_executable']         = iedb_mhc_ii_executable
            class_ii_arguments['epitope_lengths']         = [epitope_length]
            class_ii_arguments['prediction_algorithms']   = class_ii_prediction_algorithms
            class_ii_arguments['output_dir']              = per_epitope_output_dir
            class_ii_arguments['netmhc_stab']             = False
            pipeline = PvacbindPipeline(**class_ii_arguments)
            pipeline.execute()
            intermediate_output_file = os.path.join(per_epitope_output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
            output_file = os.path.join(per_epitope_output_dir, "{}.all_epitopes.final.tsv".format(args.sample_name))
            append_columns(intermediate_output_file, "{}.tsv".format(input_file), output_file)
            output_files.append(output_file)
            if epitope_length == max(args.class_ii_epitope_length):
                fasta_file = os.path.join(output_dir, "{}.fasta".format(args.sample_name))
                shutil.copy(input_file, fasta_file)
        all_epitopes_file = os.path.join(output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
        filtered_file = os.path.join(output_dir, "{}.filtered.tsv".format(args.sample_name))
        create_combined_reports(output_files, all_epitopes_file, filtered_file, True, args)
    elif len(class_ii_prediction_algorithms) == 0:
        print("No MHC class II prediction algorithms chosen. Skipping MHC class II predictions.")
    elif len(class_ii_alleles) == 0:
        print("No MHC class II alleles chosen. Skipping MHC class II predictions.")

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0 and len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        print("Creating combined reports")
        output_dir = os.path.join(base_output_dir, 'combined')
        os.makedirs(output_dir, exist_ok=True)
        file1 = os.path.join(base_output_dir, 'MHC_Class_I', "{}.all_epitopes.tsv".format(args.sample_name))
        file2 = os.path.join(base_output_dir, 'MHC_Class_II', "{}.all_epitopes.tsv".format(args.sample_name))
        combined_output_file = os.path.join(output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
        filtered_report_file = os.path.join(output_dir, "{}.filtered.tsv".format(args.sample_name))
        create_combined_reports([file1, file2], combined_output_file, filtered_report_file, False, args)

    change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
