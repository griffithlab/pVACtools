from .pvaccompare.compare_tools import *
from .pvaccompare.server import main as startserver
import argparse
import logging
from datetime import datetime
import sys
import os

logging.basicConfig(level=logging.DEBUG, format="%(message)s")


def define_parser():
    """
    Purpose:    Define arguments for the parser that the user can use
    Modifies:   Nothing
    Returns:    The parser
    """
    valid_aggregated_columns = [
        "Gene",
        "AA Change",
        "Num Passing Transcripts",
        "Best Peptide",
        "Best Transcript",
        "Num Passing Peptides",
        "IC50 MT",
        "IC50 WT",
        "%%ile MT",
        "%%ile WT",
        "RNA Expr",
        "RNA VAF",
        "DNA VAF",
        "Tier",
        "Ref Match",
    ]
    valid_unaggregated_columns = [
        "Biotype",
        "Median MT IC50 Score",
        "Median WT IC50 Score",
        "Median MT Percentile",
        "Median WT Percentile",
        "WT Epitope Seq",
        "Tumor DNA VAF",
        "Tumor RNA Depth",
        "Tumor RNA VAF",
        "Gene Expression",
        "BigMHC_EL WT Score",
        "BigMHC_EL MT Score",
        "BigMHC_IM WT Score",
        "BigMHC_IM MT Score",
        "MHCflurryEL Processing WT Score",
        "MHCflurryEL Processing MT Score",
        "MHCflurryEL Presentation WT Score",
        "MHCflurryEL Presentation MT Score",
        "MHCflurryEL Presentation WT Percentile",
        "MHCflurryEL Presentation MT Percentile",
        "MHCflurry WT IC50 Score",
        "MHCflurry MT IC50 Score",
        "MHCflurry WT Percentile",
        "MHCflurry MT Percentile",
        "MHCnuggetsI WT IC50 Score",
        "MHCnuggetsI MT IC50 Score",
        "MHCnuggetsI WT Percentile",
        "MHCnuggetsI MT Percentile",
        "NetMHC WT IC50 Score",
        "NetMHC MT IC50 Score",
        "NetMHC WT Percentile",
        "NetMHC MT Percentile",
        "NetMHCcons WT IC50 Score",
        "NetMHCcons MT IC50 Score",
        "NetMHCcons WT Percentile",
        "NetMHCcons MT Percentile",
        "NetMHCpan WT IC50 Score",
        "NetMHCpan MT IC50 Score",
        "NetMHCpan WT Percentile",
        "NetMHCpan MT Percentile",
        "NetMHCpanEL WT Score",
        "NetMHCpanEL MT Score",
        "NetMHCpanEL WT Percentile",
        "NetMHCpanEL MT Percentile",
        "PickPocket WT IC50 Score",
        "PickPocket MT IC50 Score",
        "PickPocket WT Percentile",
        "PickPocket MT Percentile",
        "SMM WT IC50 Score",
        "SMM MT IC50 Score",
        "SMM WT Percentile",
        "SMM MT Percentile",
        "SMMPMBEC WT IC50 Score",
        "SMMPMBEC MT IC50 Score",
        "SMMPMBEC WT Percentile",
        "SMMPMBEC MT Percentile",
        "DeepImmuno WT Score",
        "DeepImmuno MT Score",
        "Problematic Positions",
        "Gene of Interest",
    ]
    valid_reference_match_columns = [
        "Peptide",
        "Hit Definition",
        "Match Window",
        "Match Sequence",
    ]

    default_aggregated_columns = [
        "Num Passing Transcripts",
        "Best Peptide",
        "Best Transcript",
        "Num Passing Peptides",
        "Tier",
    ]
    default_unaggregated_columns = [
        "Biotype",
        "Median MT IC50 Score",
        "Median WT IC50 Score",
        "Median MT Percentile",
        "Median WT Percentile",
        "WT Epitope Seq",
        "Tumor DNA VAF",
        "Tumor RNA Depth",
        "Tumor RNA VAF",
        "Gene Expression",
    ]
    default_reference_match_columns = ["Peptide", "Match Window"]

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "results_folder1", nargs="?", help="Path to first results input folder"
    )
    parser.add_argument(
        "results_folder2", nargs="?", help="Path to second results input folder"
    )
    parser.add_argument(
        "--output-dir",
        default="compare_output",
        help="Specify where the output directory should be generated",
    )
    parser.add_argument(
        "--mhc-class",
        choices=["1", "2"],
        help="Specify the MHC class to run: '1' for Class I or '2' for Class II. If not specified, both classes will be processed.",
    )
    parser.add_argument(
        "--no-server",
        action="store_true",
        help="If specified, will not start the report server after the comparisons finish",
    )
    parser.add_argument(
        "--start-server",
        action="store_true",
        help="If specified, will only start the report server and will not run a comparison",
    )
    parser.add_argument(
        "--aggregated-columns",
        type=lambda s: [a for a in s.split(",")],
        default=default_aggregated_columns,
        help=f"Comma-separated columns to include in the aggregated TSV comparison, choices: {', '.join(valid_aggregated_columns)}",
    )
    parser.add_argument(
        "--unaggregated-columns",
        type=lambda s: [a for a in s.split(",")],
        default=default_unaggregated_columns,
        help=f"Comma-separated columns to include in the unaggregated TSV comparison, choices: {', '.join(valid_unaggregated_columns)}",
    )
    parser.add_argument(
        "--reference-match-columns",
        type=lambda s: [a for a in s.split(",")],
        default=default_reference_match_columns,
        help=f"Comma-separated columns to include in the reference match TSV comparison, choices: {', '.join(valid_reference_match_columns)}",
    )

    return parser


def prepare_results_folder(classes, base_output_dir):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    unique_output_dir = os.path.join(base_output_dir, f"results_{timestamp}")

    os.makedirs(unique_output_dir, exist_ok=True)

    if "1" in classes:
        os.makedirs(os.path.join(unique_output_dir, "mhc_class_i"), exist_ok=True)
    if "2" in classes:
        os.makedirs(os.path.join(unique_output_dir, "mhc_class_ii"), exist_ok=True)

    return unique_output_dir


def main(args_input=sys.argv[1:]):
    """
    Purpose:    Control function for the whole tool, calls run_comparison which calls all of the comparisons
    Modifies:   Nothing
    Returns:    None
    """
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.start_server:
        startserver()
        return

    if not args.results_folder1 or not args.results_folder2:
        logging.error(
            "Error: results_folder1 and results_folder2 are required paramaters\n"
        )
        parser.print_help()
        return

    validators.validate_aggregated_columns(args.aggregated_columns, parser)
    validators.validate_unaggregated_columns(args.unaggregated_columns, parser)
    validators.validate_reference_match_columns(args.reference_match_columns, parser)

    classes = [args.mhc_class] if args.mhc_class else ["1", "2"]
    output_dir = prepare_results_folder(classes, args.output_dir)

    for class_type in classes:
        comparison_router.run_comparison(
            class_type,
            args.results_folder1,
            args.results_folder2,
            output_dir,
            args.aggregated_columns,
            args.unaggregated_columns,
            args.reference_match_columns,
        )

    if not args.no_server:
        startserver()


if __name__ == "__main__":
    main()
