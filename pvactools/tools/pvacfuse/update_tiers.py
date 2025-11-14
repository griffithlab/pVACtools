import sys
import argparse
import tempfile

from pvactools.lib.update_tiers import UpdateTiers, PvacfuseUpdateTiers

def define_parser():
    return UpdateTiers.parser('pvacfuse')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacfuseUpdateTiers(
        args.input_file,
        binding_threshold=args.binding_threshold,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
        binding_percentile_threshold=args.binding_percentile_threshold,
        immunogenicity_percentile_threshold=args.immunogenicity_percentile_threshold,
        presentation_percentile_threshold=args.presentation_percentile_threshold,
        percentile_threshold_strategy=args.percentile_threshold_strategy,
        read_support=args.read_support,
        expn_val=args.expn_val,
        top_score_metric2=args.top_score_metric2,
    ).execute()

if __name__ == "__main__":
    main()
