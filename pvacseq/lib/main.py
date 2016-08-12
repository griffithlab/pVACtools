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
import shutil

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser("pvacseq run")

    parser.add_argument(
        "input_file",
        help="Input VCF with VEP annotations (please provide complete path)"
    )
    parser.add_argument(
        "sample_name",
        help="Name of Sample; will be used as prefix for output files"
    )
    parser.add_argument(
        "allele", type=lambda s:[a for a in s.split(',')],
        help="Allele name to predict epitope prediction. "
             + "Multiple alleles can be specified using a comma-separated list. "
             + "For a list of available alleles, use: pvacseq valid_alleles",
    )
    parser.add_argument(
        "prediction_algorithms",
        choices=PredictionClass.prediction_methods(),
        nargs="+",
        help="The epitope prediction algorithms to use",
    )
    parser.add_argument(
        "output_dir",
        help="Output directory for writing all result files"
    )
    parser.add_argument(
        "-e", "--epitope-length", type=lambda s:[int(epl) for epl in s.split(',')],
        help="Length of subpeptides(epitopes) to predict. "
             + "Multiple lengths can be specified using a comma-separated list. "
             + "Typical epitope lengths vary between 8-11. " 
             + "Required for Class II prediction algorithms",
    )
    parser.add_argument(
        "-l", "--peptide-sequence-length", type=int,
        default=21,
        help="length of the peptide sequences in the input FASTA file. Default: 21",
    )
    parser.add_argument(
        '-g', '--gene-expn-file',
        help='genes.fpkm_tracking file from Cufflinks'
    )
    parser.add_argument(
        '-i', '--transcript-expn-file',
        help='isoforms.fpkm_tracking file from Cufflinks'
    )
    parser.add_argument(
        '--net-chop-method',
        choices=lib.net_chop.methods,
        default=None,
        help="NetChop prediction method to use (\"cterm\" for C term 3.0, \"20s\" for 20S 3.0). "
             + "Default: \"cterm\" (C term 3.0)",
    )
    parser.add_argument(
        '--netmhc-stab',
        action='store_true',
        help="Run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes"
    )
    parser.add_argument(
        '-t', '--top-result-per-mutation',
        action='store_true',
        help='Output top scoring candidate per allele-length per mutation. Default: False'
    )
    parser.add_argument(
        '-m', '--top-score-metric',
        choices=['lowest', 'median'],
        default='median',
        help="The ic50 scoring metric to use when filtering epitopes. "
             + "lowest: Best MT Score - lowest MT ic50 binding score of all chosen prediction methods. "
             + "median: Median MT Score All Methods - median MT ic50 binding score of all chosen prediction methods. "
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
             + "(requiring that binding is better to the MT than WT)",
    )
    parser.add_argument(
        '--expn-val', type=int,
        default=1,
        help="Gene Expression (FPKM) Cutoff. Default: 1",
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
             + "For some resource-intensive prediction algorithms like Pickpocket and NetMHC it might be helpful to reduce this number. "
             + "Needs to be an even number.",
    )
    parser.add_argument(
        "-k", "--keep-tmp-files",
        action='store_true',
        help="Keep intermediate output files.",
    )

    args = parser.parse_args(args_input)

    PredictionClass.check_alleles_valid(args.allele)

    if "." in args.sample_name:
        sys.exit("Sample name cannot contain '.'")

    if args.fasta_size%2 != 0:
        sys.exit("The fasta size needs to be an even number")

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

    arguments = {
        'input_file'              : args.input_file,
        'sample_name'             : args.sample_name,
        'alleles'                 : args.allele,
        'gene_expn_file'          : args.gene_expn_file,
        'transcript_expn_file'    : args.transcript_expn_file,
        'net_chop_method'         : args.net_chop_method,
        'net_chop_threshold'      : args.net_chop_threshold,
        'netmhc_stab'             : args.netmhc_stab,
        'top_result_per_mutation' : args.top_result_per_mutation,
        'top_score_metric'        : args.top_score_metric,
        'binding_threshold'       : args.binding_threshold,
        'minimum_fold_change'     : args.minimum_fold_change,
        'expn_val'                : args.expn_val,
        'fasta_size'              : args.fasta_size,
        'keep_tmp_files'          : args.keep_tmp_files,
    }

    if len(class_i_prediction_algorithms) > 0:
        if args.epitope_length is None:
            sys.exit("Epitope length is required for class I binding predictions")

        print("Executing MHC Class I predictions")

        output_dir = os.path.join(base_output_dir, 'class_i')
        os.makedirs(output_dir)

        arguments['peptide_sequence_length'] = args.peptide_sequence_length
        arguments['epitope_lengths']         = args.epitope_length
        arguments['prediction_algorithms']   = class_i_prediction_algorithms
        arguments['output_dir']              = output_dir
        pipeline = MHCIPipeline(**arguments)
        pipeline.execute()

if __name__ == '__main__':
    main()
