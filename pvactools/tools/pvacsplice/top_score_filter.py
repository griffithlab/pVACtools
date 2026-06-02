import sys

from pvactools.lib.top_score_filter import PvacspliceTopScoreFilter, TopScoreFilter

def define_parser():
    return TopScoreFilter.parser('pvacsplice')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacspliceTopScoreFilter(
        args.input_file,
        args.output_file,
        top_score_metric=args.top_score_metric,
        top_score_metric2=args.top_score_metric2,
        maximum_transcript_support_level=args.maximum_transcript_support_level,
        binding_threshold=args.binding_threshold,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
        allele_specific_anchors=args.allele_specific_anchors,
        anchor_contribution_threshold=args.anchor_contribution_threshold,
        transcript_prioritization_strategy=args.transcript_prioritization_strategy,
    ).execute()

if __name__ == "__main__":
    main()
