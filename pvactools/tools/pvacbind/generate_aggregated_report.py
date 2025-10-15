import sys
import argparse
import tempfile

from pvactools.lib.aggregate_all_epitopes import PvacbindAggregateAllEpitopes
from pvactools.lib.run_utils import float_range

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
        '--binding-percentile-threshold', type=float_range(0.0,100.0),
        default=2.0,
        help="Tier epitopes in the \"Pass\" tier when the mutant allele "
             + "has a binding percentile below this value.",
    )
    parser.add_argument(
        '--immunogenicity-percentile-threshold', type=float_range(0.0,100.0),
        default=2.0,
        help="Tier epitopes in the \"Pass\" tier when the mutant allele "
             + "has a immunogenicity percentile below this value.",
    )
    parser.add_argument(
        '--presentation-percentile-threshold', type=float_range(0.0,100.0),
        default=2.0,
        help="Tier epitopes in the \"Pass\" tier when the mutant allele "
             + "has a presentation percentile below this value.",
    )
    parser.add_argument(
        '--percentile-threshold-strategy',
        choices=['conservative', 'exploratory'],
        help="Specify the candidate inclusion strategy. The 'conservative' option requires a candidate to pass the binding threshold and all percentile thresholds (default)."
             + " The 'exploratory' option requires a candidate to pass at the binding threshold or one of the percentile thresholds.",
        default="conservative",
    )
    parser.add_argument(
        '--aggregate-inclusion-binding-threshold', type=int,
        help="Threshold for including epitopes when creating the aggregate report",
        default=5000,
    )
    parser.add_argument(
        '--aggregate-inclusion-count-limit', type=int,
        help="Limit neoantigen candidates included in the aggregate report to only the best n candidates per variant.",
        default=15,
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
        '-m2', '--top-score-metric2',
        choices=['ic50','percentile'],
        default='ic50',
        help="Whether to use median/best IC50 or to use median/best percentile score when determining the best peptide in the aggregated report. "
             + "This parameter is also used to influence the primary sorting criteria in the aggregated report for the candidates within each tier."
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    tmp_fh = tempfile.NamedTemporaryFile()

    print("Creating Aggregated Report")
    PvacbindAggregateAllEpitopes(
        args.input_file,
        args.output_file,
        binding_threshold=args.binding_threshold,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
        binding_percentile_threshold=args.binding_percentile_threshold,
        immunogenicity_percentile_threshold=args.immunogenicity_percentile_threshold,
        presentation_percentile_threshold=args.presentation_percentile_threshold,
        percentile_threshold_strategy=args.percentile_threshold_strategy,
        top_score_metric=args.top_score_metric,
        aggregate_inclusion_binding_threshold=args.aggregate_inclusion_binding_threshold,
        aggregate_inclusion_count_limit=args.aggregate_inclusion_count_limit,
    ).execute()
    print("Completed")

if __name__ == '__main__':
    main()
