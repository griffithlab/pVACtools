import sys
import argparse
import os
from lib.prediction_class import *
from lib.pipeline import *
from tools.pvacseq.config_files import additional_input_file_list_options
from lib.run_argument_parser import *
from lib.condense_final_report import *
from lib.rank_epitopes import *
from lib.binding_filter import *
from lib.top_score_filter import *

def define_parser():
    return PvacfuseRunArgumentParser().parser

def combine_reports(input_files, output_file):
    write_headers = True
    with open(output_file, 'w') as fout:
        writer = csv.writer(fout)
        for filename in input_files:
            with open(filename) as fin:
                reader = csv.reader(fin)
                headers = next(reader)
                if write_headers:
                    write_headers = False  # Only write headers once.
                    writer.writerow(headers)
                writer.writerows(reader)  # Write all remaining rows.

def binding_filter(input_file, output_dir, args):
    output_file = os.path.join(output_dir, "{}.filtered.binding.tsv".format(args.sample_name))
    print("Running Binding Filters")
    BindingFilter(
        input_file,
        output_file,
        args.binding_threshold,
        0,
        args.top_score_metric,
        args.exclude_NAs,
    ).execute()
    print("Completed")
    return output_file

def top_result_filter(coverage_filter_output_file, output_dir, args):
    output_file = os.path.join(output_dir, "{}.filtered.top.tsv".format(args.sample_name))
    print("Running Top Score Filter")
    TopScoreFilter(coverage_filter_output_file, output_file, args.top_score_metric).execute()
    print("Completed")
    return output_file

def condensed_report(final_output_file, output_dir, args):
    output_file = os.path.join(output_dir, "{}.final.condensed.tsv".format(args.sample_name))
    print("Creating condensed final report")
    CondenseFinalReport(final_output_file, output_file, args.top_score_metric).execute()
    print("Completed")
    return output_file

def rank_epitopes(condensed_report_output_file, output_dir, args):
    output_file = os.path.join(output_dir, "{}.filtered.condensed.ranked.tsv".format(args.sample_name))
    print("Ranking neoepitopes")
    RankEpitopes(condensed_report_output_file, output_file).execute()
    print("Completed")
    return output_file

def create_combined_reports(base_output_dir, args):
    output_dir = os.path.join(base_output_dir, 'combined')
    os.makedirs(output_dir, exist_ok=True)

    file1 = os.path.join(base_output_dir, 'MHC_Class_I', "{}.all_epitopes.tsv".format(args.sample_name))
    file2 = os.path.join(base_output_dir, 'MHC_Class_II', "{}.all_epitopes.tsv".format(args.sample_name))
    combined_output_file = os.path.join(output_dir, "{}.all_epitopes.tsv".format(args.sample_name))
    combine_reports([file1, file2], combined_output_file)

    binding_filter_output_file = binding_filter(combined_output_file, output_dir, args)
    top_result_filter_output_file = top_result_filter(binding_filter_output_file, output_dir, args)
    final_output_file = os.path.join(output_dir, "{}.filtered.tsv".format(args.sample_name))
    shutil.copy(top_result_filter_output_file, final_output_file)
    condensed_report_output_file = condensed_report(final_output_file, output_dir, args)
    ranked_output_file = rank_epitopes(condensed_report_output_file, output_dir, args)
    for file_name in [
        binding_filter_output_file,
        top_result_filter_output_file,
        condensed_report_output_file,
    ]:
        os.unlink(file_name)
    print("\nDone: Pipeline finished successfully. File {} contains ranked list of filtered putative neoantigens for class I and class II predictions.\n".format(ranked_output_file))

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if "." in args.sample_name:
        sys.exit("Sample name cannot contain '.'")

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

    input_file_type = 'bedpe'
    base_output_dir = os.path.abspath(args.output_dir)

    class_i_prediction_algorithms = []
    class_ii_prediction_algorithms = []
    for prediction_algorithm in sorted(args.prediction_algorithms):
        prediction_class = globals()[prediction_algorithm]
        prediction_class_object = prediction_class()
        if isinstance(prediction_class_object, MHCI):
            class_i_prediction_algorithms.append(prediction_algorithm)
        elif isinstance(prediction_class_object, MHCII):
            class_ii_prediction_algorithms.append(prediction_algorithm)

    class_i_alleles = []
    class_ii_alleles = []
    for allele in sorted(set(args.allele)):
        valid = 0
        if allele in MHCI.all_valid_allele_names():
            class_i_alleles.append(allele)
            valid = 1
        if allele in MHCII.all_valid_allele_names():
            class_ii_alleles.append(allele)
            valid = 1
        if not valid:
            print("Allele %s not valid. Skipping." % allele)

    shared_arguments = {
        'input_file'                : args.input_file,
        'input_file_type'           : input_file_type,
        'sample_name'               : args.sample_name,
        'top_score_metric'          : args.top_score_metric,
        'binding_threshold'         : args.binding_threshold,
        'net_chop_method'           : args.net_chop_method,
        'net_chop_threshold'        : args.net_chop_threshold,
        'additional_report_columns' : args.additional_report_columns,
        'fasta_size'                : args.fasta_size,
        'iedb_retries'              : args.iedb_retries,
        'downstream_sequence_length': downstream_sequence_length,
        'keep_tmp_files'            : args.keep_tmp_files,
    }

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0:
        if args.epitope_length is None:
            sys.exit("Epitope length is required for class I binding predictions")

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
        class_i_arguments['peptide_sequence_length'] = args.peptide_sequence_length
        class_i_arguments['iedb_executable']         = iedb_mhc_i_executable
        class_i_arguments['epitope_lengths']         = args.epitope_length
        class_i_arguments['prediction_algorithms']   = class_i_prediction_algorithms
        class_i_arguments['output_dir']              = output_dir
        class_i_arguments['netmhc_stab']             = args.netmhc_stab
        pipeline = MHCIPipeline(**class_i_arguments)
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
        class_ii_arguments['alleles']               = class_ii_alleles
        class_ii_arguments['prediction_algorithms'] = class_ii_prediction_algorithms
        class_ii_arguments['iedb_executable']       = iedb_mhc_ii_executable
        class_ii_arguments['output_dir']            = output_dir
        class_ii_arguments['netmhc_stab']           = False
        pipeline = MHCIIPipeline(**class_ii_arguments)
        pipeline.execute()
    elif len(class_ii_prediction_algorithms) == 0:
        print("No MHC class II prediction algorithms chosen. Skipping MHC class II predictions.")
    elif len(class_ii_alleles) == 0:
        print("No MHC class II alleles chosen. Skipping MHC class II predictions.")

    if len(class_i_prediction_algorithms) > 0 and len(class_i_alleles) > 0 and len(class_ii_prediction_algorithms) > 0 and len(class_ii_alleles) > 0:
        print("Creating combined reports")
        create_combined_reports(base_output_dir, args)

if __name__ == '__main__':
    main()
