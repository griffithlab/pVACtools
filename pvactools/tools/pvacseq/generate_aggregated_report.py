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
    parser.add_argument(
        '-b', '--binding-threshold', type=int,
        help="Tier epitopes in the \"Pass\" tier when the mutant allele "
             + "has ic50 binding scores below this value and in the \"Relaxed\" tier when the mutation allele has ic50 binding scores below double this value.",
        default=500
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
        trna_vaf=args.trna_vaf,
        trna_cov=args.trna_cov,
        expn_val=args.expn_val
    ).execute()
    print("Completed")

if __name__ == '__main__':
    main()
