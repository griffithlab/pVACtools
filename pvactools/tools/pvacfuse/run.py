import sys
import argparse
import os
import shutil

from pvactools.lib.prediction_class import *
from pvactools.lib.pipeline import PvacbindPipeline
from pvactools.lib.run_argument_parser import PvacfuseRunArgumentParser
from pvactools.lib.post_processor import PostProcessor
import pvactools.tools.pvacfuse.generate_protein_fasta
import pvactools.lib.run_utils

def define_parser():
    return PvacfuseRunArgumentParser().parser

def create_net_class_report(files, all_epitopes_output_file, filtered_report_file, args, run_params):
    for file_name in files:
        if not os.path.exists(file_name):
            print("File {} doesn't exist. Aborting.".format(file_name))
            return

    pvactools.lib.run_utils.combine_reports(files, all_epitopes_output_file)

    post_processing_params = vars(args).copy()
    post_processing_params['input_file'] = all_epitopes_output_file
    post_processing_params['filtered_report_file'] = filtered_report_file
    post_processing_params['minimum_fold_change'] = None
    post_processing_params['run_coverage_filter'] = False
    post_processing_params['run_transcript_support_level_filter'] = False
    post_processing_params['run_manufacturability_metrics'] = True
    if run_params['net_chop_method']:
        post_processing_params['run_net_chop'] = True
        post_processing_params['net_chop_fasta'] = run_params['net_chop_fasta']
    else:
        post_processing_params['run_net_chop'] = False
    post_processing_params['run_netmhc_stab'] = True if run_params['netmhc_stab'] else False
    post_processing_params['fasta'] = run_params['fasta']
    post_processing_params['species'] = run_params['species']
    post_processing_params['file_type'] = 'pVACfuse'

    PostProcessor(**post_processing_params).execute()

def create_combined_reports(files, all_epitopes_output_file, filtered_report_file, run_manufacturability_metrics, args):
    for file_name in files:
        if not os.path.exists(file_name):
            print("File {} doesn't exist. Aborting.".format(file_name))
            return

    pvactools.lib.run_utils.combine_reports(files, all_epitopes_output_file)

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

def generate_fasta(args, output_dir, epitope_length, epitope_flank_length=0, net_chop_fasta=False):
    if net_chop_fasta:
        per_epitope_output_dir = None
        output_file = os.path.join(output_dir, "{}.net_chop.fa".format(args.sample_name))
    else:
        per_epitope_output_dir = os.path.join(output_dir, str(epitope_length))
        os.makedirs(per_epitope_output_dir, exist_ok=True)
        output_file = os.path.join(per_epitope_output_dir, "{}.fa".format(args.sample_name))
    params = [
        args.input_file,
        str(epitope_flank_length + epitope_length - 1),
        output_file,
    ]
    if args.downstream_sequence_length is not None:
        params.extend(["-d", str(args.downstream_sequence_length)])
    else:
        params.extend(["-d", 'full'])
    pvactools.tools.pvacfuse.generate_protein_fasta.main(params, save_tsv_file=True)
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

    (class_i_prediction_algorithms, class_ii_prediction_algorithms) = pvactools.lib.run_utils.split_algorithms(args.prediction_algorithms)
    (class_i_alleles, class_ii_alleles, species) = pvactools.lib.run_utils.split_alleles(args.allele)
    class_ii_alleles = pvactools.lib.run_utils.combine_class_ii_alleles(class_ii_alleles)

    shared_arguments = {
        'input_file_type'           : 'fasta',
        'sample_name'               : args.sample_name,
        'top_score_metric'          : args.top_score_metric,
        'binding_threshold'         : args.binding_threshold,
        'percentile_threshold'      : args.percentile_threshold,
        'allele_specific_binding_thresholds': args.allele_specific_binding_thresholds,
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
        'blastp_path'               : args.blastp_path,
        'blastp_db'                 : args.blastp_db,
        'run_post_processor'        : False,
        'exclude_NAs'               : args.exclude_NAs,
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

    all_params = {
        'I': {
            'iedb_executable': iedb_mhc_i_executable,
            'prediction_algorithms': class_i_prediction_algorithms,
            'alleles': class_i_alleles,
            'epitope_lengths': args.class_i_epitope_length,
            'netmhc_stab': args.netmhc_stab
        },
        'II': {
            'iedb_executable': iedb_mhc_ii_executable,
            'prediction_algorithms': class_ii_prediction_algorithms,
            'alleles': class_ii_alleles,
            'epitope_lengths': args.class_ii_epitope_length,
            'netmhc_stab': False
        }
    }

    for (mhc_class, params) in all_params.items():
        prediction_algorithms = params['prediction_algorithms']
        alleles = params['alleles']
        epitope_lengths = params['epitope_lengths']
        iedb_executable = params['iedb_executable']
        netmhc_stab = params['netmhc_stab']

        if len(prediction_algorithms) > 0 and len(alleles) > 0:
            print("Executing MHC Class {} predictions".format(mhc_class))

            output_dir = os.path.join(base_output_dir, 'MHC_Class_{}'.format(mhc_class))
            os.makedirs(output_dir, exist_ok=True)

            output_files = []
            run_arguments = shared_arguments.copy()
            run_arguments['alleles']               = alleles
            run_arguments['iedb_executable']       = iedb_executable
            run_arguments['prediction_algorithms'] = prediction_algorithms
            run_arguments['netmhc_stab']           = netmhc_stab

            for epitope_length in epitope_lengths:
                (input_file, per_epitope_output_dir) = generate_fasta(args, output_dir, epitope_length)
                if os.path.getsize(input_file) == 0:
                    print("The intermediate FASTA file for epitope length {} is empty. Please check that the input AGfusion directory contains fusion entries with `*_protein.fa` files. Fusion entries without this file cannot be processed by pVACfuse.".format(epitope_length))
                    continue

                run_arguments['input_file']              = input_file
                run_arguments['epitope_lengths']         = [epitope_length]
                run_arguments['output_dir']              = per_epitope_output_dir
                pipeline = PvacbindPipeline(**run_arguments)
                pipeline.execute()
                intermediate_output_file = os.path.join(per_epitope_output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
                output_file = os.path.join(per_epitope_output_dir, "{}.all_epitopes.final.tsv".format(args.sample_name))
                append_columns(intermediate_output_file, "{}.tsv".format(input_file), output_file)
                output_files.append(output_file)
                if epitope_length == max(epitope_lengths):
                    # copy fasta to output dir
                    fasta_file = os.path.join(output_dir, "{}.fasta".format(args.sample_name))
                    shutil.copy(input_file, fasta_file)
                    run_arguments['fasta'] = fasta_file
                    # generate and copy net_chop fasta to output dir if specified
                    if args.net_chop_method:
                        epitope_flank_length = 9
                        (net_chop_fasta, _) = generate_fasta(args, output_dir, epitope_length, epitope_flank_length, net_chop_fasta=True)
                        run_arguments['net_chop_fasta'] = net_chop_fasta
            if len(output_files) > 0:
                all_epitopes_file = os.path.join(output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
                filtered_file = os.path.join(output_dir, "{}.filtered.tsv".format(args.sample_name))
                #!!! make below call to create_net_class_report
                #create_combined_reports(output_files, all_epitopes_file, filtered_file, True, args)
                create_net_class_report(output_files, all_epitopes_file, filtered_file, args, run_arguments)
        elif len(prediction_algorithms) == 0:
            print("No MHC class {} prediction algorithms chosen. Skipping MHC class I predictions.".format(mhc_class))
        elif len(alleles) == 0:
            print("No MHC class{} alleles chosen. Skipping MHC class I predictions.".format(mhc_class))

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0 and len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        print("Creating combined reports")
        output_dir = os.path.join(base_output_dir, 'combined')
        os.makedirs(output_dir, exist_ok=True)
        file1 = os.path.join(base_output_dir, 'MHC_Class_I', "{}.all_epitopes.tsv".format(args.sample_name))
        file2 = os.path.join(base_output_dir, 'MHC_Class_II', "{}.all_epitopes.tsv".format(args.sample_name))
        combined_output_file = os.path.join(output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
        filtered_report_file = os.path.join(output_dir, "{}.filtered.tsv".format(args.sample_name))
        create_combined_reports([file1, file2], combined_output_file, filtered_report_file, False, args)

    pvactools.lib.run_utils.change_permissions_recursive(base_output_dir, 0o755, 0o644)

if __name__ == '__main__':
    main()
