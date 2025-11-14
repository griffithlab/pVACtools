import sys
import argparse
import tempfile

from pvactools.lib.update_tiers import UpdateTiers, PvacseqUpdateTiers

def define_parser():
    return UpdateTiers.parser('pvacseq')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacseqUpdateTiers(
        args.input_file,
        args.vaf_clonal,
        binding_threshold=args.binding_threshold,
        trna_vaf=args.trna_vaf,
        trna_cov=args.trna_cov,
        expn_val=args.expn_val,
        transcript_prioritization_strategy=args.transcript_prioritization_strategy,
        maximum_transcript_support_level=args.maximum_transcript_support_level,
        binding_percentile_threshold=args.binding_percentile_threshold,
        immunogenicity_percentile_threshold=args.immunogenicity_percentile_threshold,
        presentation_percentile_threshold=args.presentation_percentile_threshold,
        percentile_threshold_strategy=args.percentile_threshold_strategy,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
        allele_specific_anchors=args.allele_specific_anchors,
        anchor_contribution_threshold=args.anchor_contribution_threshold,
        top_score_metric2=args.top_score_metric2,
        metrics_file=args.metrics_file,
    ).execute()

if __name__ == "__main__":
    main()
