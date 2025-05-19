import sys
import argparse
import tempfile

from pvactools.lib.update_tiers import UpdateTiers, PvacbindUpdateTiers

def define_parser():
    return UpdateTiers.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacbindUpdateTiers(
        args.input_file,
        binding_threshold=args.binding_threshold,
        percentile_threshold=args.percentile_threshold,
        percentile_threshold_strategy=args.percentile_threshold_strategy,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
    ).execute()

if __name__ == "__main__":
    main()
