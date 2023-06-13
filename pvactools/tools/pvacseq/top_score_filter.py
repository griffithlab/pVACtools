import sys

from pvactools.lib.top_score_filter import PvacseqTopScoreFilter, TopScoreFilter

def define_parser():
    return TopScoreFilter.parser('pvacseq')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacseqTopScoreFilter(
        args.input_file,
        args.output_file,
        top_score_metric=args.top_score_metric,
        binding_threshold=args.binding_threshold,
        allele_specific_binding_thresholds=args.allele_specific_binding_thresholds,
        allele_specific_anchors=args.allele_specific_anchors,
        anchor_contribution_threshold=args.anchor_contribution_threshold,
    ).execute()

if __name__ == "__main__":
    main()
