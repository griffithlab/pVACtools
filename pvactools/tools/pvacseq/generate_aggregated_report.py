import sys
import argparse
import tempfile

from pvactools.lib.aggregate_all_epitopes import PvacseqAggregateAllEpitopes

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
    PvacseqAggregateAllEpitopes(args.input_file, args.output_file, args.tumor_purity).execute()
    print("Completed")

if __name__ == '__main__':
    main()
