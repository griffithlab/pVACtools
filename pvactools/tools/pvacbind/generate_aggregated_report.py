import sys
import argparse
import tempfile

from pvactools.lib.aggregate_all_epitopes import PvacbindAggregateAllEpitopes

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacbind generate_aggregated_report",
        description="Generate an aggregated report from a pVACbind .all_epitopes.tsv report file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input_file",
        help="A pVACbind .all_epitopes.tsv report file"
    )
    parser.add_argument(
        "output_file",
        help="The file path to write the aggregated report tsv to"
    )
    parser.add_argument(
        '-b', '--binding-threshold', type=int,
        help="Tier epitopes in the \"Pass\" tier when the mutant allele "
             + "has ic50 binding scores below this value and in the \"Relaxed\" tier when the mutant allele has ic50 binding scores below double this value.",
        default=500
    )
    parser.add_argument(
        '--allele-specific-binding-thresholds',
        help="Use allele-specific binding thresholds. To print the allele-specific binding thresholds run `pvacbind allele_specific_cutoffs`. "
             + "If an allele does not have a special threshold value, the `--binding-threshold` value will be used.",
        default=False,
        action='store_true',
    )
    parser.add_argument(
        '--percentile-threshold', type=float,
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

    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    tmp_fh = tempfile.NamedTemporaryFile()

    print("Creating Aggreggated Report")
    PvacbindAggregateAllEpitopes(
        args.input_file,
        args.output_file,
        binding_threshold=args.binding_threshold,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
        percentile_threshold=args.percentile_threshold,
        top_score_metric=args.top_score_metric,
        aggregate_inclusion_binding_threshold=args.aggregate_inclusion_binding_threshold,
    ).execute()
    print("Completed")

if __name__ == '__main__':
    main()
