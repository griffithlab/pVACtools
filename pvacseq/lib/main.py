import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import argparse
import os

try:
    from .. import lib
except ValueError:
    import lib
from lib.prediction_class import *
from lib.pipeline import *
from lib.config_files import additional_input_file_list_options

import shutil
import yaml

def parse_additional_input_file_list(additional_input_file_list):
    if additional_input_file_list:
        with open(additional_input_file_list, 'r') as stream:
            additional_input_files = yaml.load(stream)
        for additional_input_file in additional_input_files:
            if additional_input_file not in additional_input_file_list_options().keys():
                sys.exit("%s not a valid key in the additional_input_file_list" % additional_input_file)
    else:
        additional_input_files = {}
    for additional_input_file_list_option in additional_input_file_list_options().keys():
        if additional_input_file_list_option not in additional_input_files:
            additional_input_files[additional_input_file_list_option] = None
    return additional_input_files

def define_parser():
    parser = argparse.ArgumentParser("pvacseq run")

    parser.add_argument(
        "input_file",
        help="A VEP-annotated single-sample VCF containing transcript, Wildtype protein sequence, and Downstream protein sequence information"
    )
    parser.add_argument(
        "sample_name",
        help="The name of the sample being processed. This will be used as a prefix for output files"
    )
    parser.add_argument(
        "allele", type=lambda s:[a for a in s.split(',')],
        help="Name of the allele to use for epitope prediction. "
             + "Multiple alleles can be specified using a comma-separated list. "
             + "For a list of available alleles, use: `pvacseq valid_alleles`",
    )
    parser.add_argument(
        "prediction_algorithms",
        choices=PredictionClass.prediction_methods(),
        nargs="+",
        help="The epitope prediction algorithms to use. Multiple prediction algorithms can be specified, separated by spaces",
    )
    parser.add_argument(
        "output_dir",
        help="The directory for writing all result files"
    )
    parser.add_argument(
        "-e", "--epitope-length", type=lambda s:[int(epl) for epl in s.split(',')],
        help="Length of subpeptides (neoepitopes) to predict. "
             + "Multiple epitope lengths can be specified using a comma-separated list. "
             + "Typical epitope lengths vary between 8-11. " 
             + "Required for Class I prediction algorithms",
    )
    parser.add_argument(
        "-l", "--peptide-sequence-length", type=int,
        default=21,
        help="Length of the peptide sequence to use when creating the FASTA. Default: 21",
    )
    parser.add_argument(
        "-i", "--additional-input-file-list",
        help="yaml file of additional files to be used as inputs, e.g. cufflinks output files. "
             + "For an example of this yaml file run `pvacseq config_files additional_input_file_list`."
    )
    parser.add_argument(
        '--net-chop-method',
        choices=lib.net_chop.methods,
        default=None,
        help="NetChop prediction method to use (\"cterm\" for C term 3.0, \"20s\" for 20S 3.0). ",
    )
    parser.add_argument(
        '--netmhc-stab',
        action='store_true',
        help="Run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes"
    )
    parser.add_argument(
        '-t', '--top-result-per-mutation',
        action='store_true',
        help='Output only the top scoring result for each allele-peptide length combination for each variant. Default: False'
    )
    parser.add_argument(
        '-m', '--top-score-metric',
        choices=['lowest', 'median'],
        default='median',
        help="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. "
             + "lowest: Best MT Score/Corresponding Fold Change - lowest MT ic50 binding score/corresponding fold change of all chosen prediction methods. "
             + "median: Median MT Score/Median Fold Change - median MT ic50 binding score/fold change of all chosen prediction methods. "
             + "Default: median"
        )
    parser.add_argument(
        "-b","--binding-threshold", type=int,
        default=500,
        help="Report only epitopes where the mutant allele has ic50 binding scores below this value. Default: 500",
    )
    parser.add_argument(
        "-c", "--minimum-fold-change", type=int,
        default=0,
        help="Minimum fold change between mutant binding score and wild-type score. "
             + "The default is 0, which filters no results, but 1 is often a sensible choice "
             + "(requiring that binding is better to the MT than WT). Default: 0",
    )
    parser.add_argument(
        '--normal-cov', type=int,
        help="Normal Coverage Cutoff. Sites above this cutoff will be considered. " +
        "Default: 5",
        default=5
    )
    parser.add_argument(
        '--tdna-cov', type=int,
        help="Tumor DNA Coverage Cutoff. Sites above this cutoff will be considered. " +
        "Default: 10",
        default=10
    )
    parser.add_argument(
        '--trna-cov', type=int,
        help="Tumor RNA Coverage Cutoff. Sites above this cutoff will be considered. " +
        "Default: 10",
        default=10
    )
    parser.add_argument(
        '--normal-vaf', type=int,
        help="Normal VAF Cutoff. Sites BELOW this cutoff in normal will be considered. " +
        "Default: 2",
        default=2
    )
    parser.add_argument(
        '--tdna-vaf', type=int,
        help="Tumor DNA VAF Cutoff. Sites above this cutoff will be considered. " +
        "Default: 40",
        default=40
    )
    parser.add_argument(
        '--trna-vaf', type=int,
        help="Tumor RNA VAF Cutoff. Sites above this cutoff will be considered. " +
        "Default: 40",
        default=40
    )
    parser.add_argument(
        '--expn-val', type=int,
        default=1,
        help="Gene and Transcript Expression cutoff. Sites above this cutoff will be considered. Default: 1",
    )
    parser.add_argument(
        '--net-chop-threshold', type=float,
        default=0.5,
        help="NetChop prediction threshold. Default: 0.5",
    )
    parser.add_argument(
        "-s", "--fasta-size",type=int,
        default=200,
        help="Number of fasta entries per IEDB request. "
             + "For some resource-intensive prediction algorithms like Pickpocket and NetMHCpan it might be helpful to reduce this number. "
             + "Needs to be an even number.",
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default='1000',
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence. Default: 1000"
    )
    parser.add_argument(
        "-k", "--keep-tmp-files",
        action='store_true',
        help="Keep intermediate output files. This migt be useful for debugging purposes.",
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PredictionClass.check_alleles_valid(args.allele)

    if "." in args.sample_name:
        sys.exit("Sample name cannot contain '.'")

    if args.fasta_size%2 != 0:
        sys.exit("The fasta size needs to be an even number")

    if args.downstream_sequence_length == 'full':
        downstream_sequence_length = None
    elif args.downstream_sequence_length.isdigit():
        downstream_sequence_length = args.downstream_sequence_length
    else:
        sys.exit("The downstream sequence length needs to be a positive integer or 'full'")

    base_output_dir = os.path.abspath(args.output_dir)

    class_i_prediction_algorithms = []
    class_ii_prediction_algorithms = []
    for prediction_algorithm in args.prediction_algorithms:
        prediction_class = globals()[prediction_algorithm]
        prediction_class_object = prediction_class()
        if isinstance(prediction_class_object, MHCI):
            class_i_prediction_algorithms.append(prediction_algorithm)
        elif isinstance(prediction_class_object, MHCII):
            class_ii_prediction_algorithms.append(prediction_algorithm)

    shared_arguments = {
        'input_file'                : args.input_file,
        'sample_name'               : args.sample_name,
        'alleles'                   : args.allele,
        'top_result_per_mutation'   : args.top_result_per_mutation,
        'top_score_metric'          : args.top_score_metric,
        'binding_threshold'         : args.binding_threshold,
        'minimum_fold_change'       : args.minimum_fold_change,
        'net_chop_method'           : args.net_chop_method,
        'net_chop_threshold'        : args.net_chop_threshold,
        'normal_cov'                : args.normal_cov,
        'normal_vaf'                : args.normal_vaf,
        'tdna_cov'                  : args.tdna_cov,
        'tdna_vaf'                  : args.tdna_vaf,
        'trna_cov'                  : args.trna_cov,
        'trna_vaf'                  : args.trna_vaf,
        'expn_val'                  : args.expn_val,
        'fasta_size'                : args.fasta_size,
        'downstream_sequence_length': downstream_sequence_length,
        'keep_tmp_files'            : args.keep_tmp_files,
    }
    additional_input_files = parse_additional_input_file_list(args.additional_input_file_list)
    shared_arguments.update(additional_input_files)

    if len(class_i_prediction_algorithms) > 0:
        if args.epitope_length is None:
            sys.exit("Epitope length is required for class I binding predictions")

        print("Executing MHC Class I predictions")

        output_dir = os.path.join(base_output_dir, 'MHC_Class_I')
        os.makedirs(output_dir, exist_ok=True)

        class_i_arguments = shared_arguments.copy()
        class_i_arguments['peptide_sequence_length'] = args.peptide_sequence_length
        class_i_arguments['epitope_lengths']         = args.epitope_length
        class_i_arguments['prediction_algorithms']   = class_i_prediction_algorithms
        class_i_arguments['output_dir']              = output_dir
        class_i_arguments['netmhc_stab']             = args.netmhc_stab
        pipeline = MHCIPipeline(**class_i_arguments)
        pipeline.execute()

    if len(class_ii_prediction_algorithms) > 0:
        print("Executing MHC Class II predictions")

        output_dir = os.path.join(base_output_dir, 'MHC_Class_II')
        os.makedirs(output_dir, exist_ok=True)

        class_ii_arguments = shared_arguments.copy()
        class_ii_arguments['prediction_algorithms'] = class_ii_prediction_algorithms
        class_ii_arguments['output_dir']            = output_dir
        class_ii_arguments['netmhc_stab']           = False
        pipeline = MHCIIPipeline(**class_ii_arguments)
        pipeline.execute()

if __name__ == '__main__':
    main()
