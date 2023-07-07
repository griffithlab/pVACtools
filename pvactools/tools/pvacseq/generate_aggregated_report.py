import sys
import argparse
import tempfile

from pvactools.lib.aggregate_all_epitopes import PvacseqAggregateAllEpitopes
from pvactools.lib.run_utils import *

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacseq generate_aggregated_report",
        description="Generate an aggregated report from a pVACseq .all_epitopes.tsv report file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input_file",
        help="A pVACseq .all_epitopes.tsv report file"
    )
    parser.add_argument(
        "output_file",
        help="The file path to write the aggregated report tsv to"
    )
    parser.add_argument(
        "--tumor-purity",
        help="Value between 0 and 1 indicating the fraction of tumor cells in the tumor sample. Information is used during aggregate report creation for a simple estimation of whether variants are subclonal or clonal based on VAF. If not provided, purity is estimated directly from the VAFs.",
        type=float,
    )
    parser.add_argument(
        '-b', '--binding-threshold', type=int,
        help="Tier epitopes in the \"Pass\" tier when the mutant allele "
             + "has ic50 binding scores below this value.",
        default=500
    )
    parser.add_argument(
        '--allele-specific-binding-thresholds',
        help="Use allele-specific binding thresholds. To print the allele-specific binding thresholds run `pvacseq allele_specific_cutoffs`. "
             + "If an allele does not have a special threshold value, the `--binding-threshold` value will be used.",
        default=False,
        action='store_true',
    )
    parser.add_argument(
        '--percentile-threshold', type=float_range(0.0,100.0),
        help="When set, tier epitopes in the \"Pass\" tier when the mutant allele "
             + "has percentile scores below this value and in the \"Relaxed\" tier "
             + "when the mutant allele has percentile scores below double this value.",
    )
    parser.add_argument(
        '--aggregate-inclusion-binding-threshold', type=int,
        help="Threshold for including epitopes when creating the aggregate report",
        default=5000,
    )
    parser.add_argument(
        '-m', '--top-score-metric',
        choices=['lowest', 'median'],
        default='median',
        help="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. "
             + "lowest: Use the best MT Score and Corresponding Fold Change (i.e. the lowest MT ic50 binding score and corresponding fold change of all chosen prediction methods). "
             + "median: Use the median MT Score and Median Fold Change (i.e. the  median MT ic50 binding score and fold change of all chosen prediction methods)."
    )
    parser.add_argument(
        '--trna-vaf', type=float,
        help="Tumor RNA VAF Cutoff. Used to calculate the allele expression cutoff for tiering.",
        default=0.25
    )
    parser.add_argument(
        '--trna-cov', type=int,
        help="Tumor RNA Coverage Cutoff. Used as a cutoff for tiering.",
        default=10
    )
    parser.add_argument(
        '--expn-val', type=float,
        default=1.0,
        help="Gene and Expression cutoff. Used to calculate the allele expression cutoff for tiering.",
    )
    parser.add_argument(
        "--maximum-transcript-support-level", type=int,
        help="The threshold to use for filtering epitopes on the Ensembl transcript support level (TSL). "
        +"Transcript support level needs to be <= this cutoff to be included in most tiers.",
        default=1,
        choices=[1,2,3,4,5]
    )
    parser.add_argument(
        "--allele-specific-anchors",
        help="Use allele-specific anchor positions when tiering epitopes in the aggregate report. This option "
             + "is available for 8, 9, 10, and 11mers and only for HLA-A, B, and C alleles. If this option is "
             + "not enabled or as a fallback for unsupported lengths and alleles, the default positions of 1, "
             + "2, epitope length - 1, and epitope length are used. Please see https://doi.org/10.1101/2020.12.08.416271 "
             + "for more details.",
        default=False,
        action='store_true',
    )
    parser.add_argument(
        "--anchor-contribution-threshold", type=float_range(0.5,0.9),
        help="For determining allele-specific anchors, each position is assigned a score based on how binding is "
             + "influenced by mutations. From these scores, the relative contribution of each position to the "
             + "overall binding is calculated. Starting with the highest relative contribution, positions whose "
             + "scores together account for the selected contribution threshold are assigned as anchor locations. "
             + " As a result, a higher threshold leads to the inclusion of more positions to be considered anchors.",
        default=0.8
    )

    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.tumor_purity is not None:
        if args.tumor_purity > 1:
            raise Exception("--tumor-purity must be a float between 0 and 1. Value too large: {}".format(args.tumor_purity))
        elif args.tumor_purity < 0:
            raise Exception("--tumor-purity must be a float between 0 and 1. Value too small: {}".format(args.tumor_purity))

    tmp_fh = tempfile.NamedTemporaryFile()

    print("Creating Aggreggated Report")
    PvacseqAggregateAllEpitopes(
        args.input_file,
        args.output_file,
        tumor_purity=args.tumor_purity,
        binding_threshold=args.binding_threshold,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
        percentile_threshold=args.percentile_threshold,
        trna_vaf=args.trna_vaf,
        trna_cov=args.trna_cov,
        expn_val=args.expn_val,
        maximum_transcript_support_level=args.maximum_transcript_support_level,
        top_score_metric=args.top_score_metric,
        allele_specific_anchors=args.allele_specific_anchors,
        anchor_contribution_threshold=args.anchor_contribution_threshold,
        aggregate_inclusion_binding_threshold=args.aggregate_inclusion_binding_threshold,
    ).execute()
    print("Completed")

if __name__ == '__main__':
    main()
