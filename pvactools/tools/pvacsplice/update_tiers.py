import sys
import argparse
import tempfile

from pvactools.lib.update_tiers import UpdateTiers, PvacspliceUpdateTiers

def define_parser():
    return UpdateTiers.parser('pvacsplice')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacspliceUpdateTiers(
        args.input_file,
        args.vaf_clonal,
        binding_threshold=args.binding_threshold,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
        percentile_threshold=args.percentile_threshold,
        percentile_threshold_strategy=args.percentile_threshold_strategy,
        trna_vaf=args.trna_vaf,
        trna_cov=args.trna_cov,
        expn_val=args.expn_val,
        transcript_prioritization_strategy=args.transcript_prioritization_strategy,
        maximum_transcript_support_level=args.maximum_transcript_support_level,
        top_score_metric2=args.top_score_metric2,
    ).execute()

if __name__ == "__main__":
    main()
