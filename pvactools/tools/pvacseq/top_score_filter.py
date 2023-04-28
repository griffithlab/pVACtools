import sys

from pvactools.lib.top_score_filter import TopScoreFilter

def define_parser():
    return TopScoreFilter.parser('pvacseq')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    TopScoreFilter(
        args.input_file,
        args.output_file,
        args.top_score_metric,
        args.binding_threshold,
        args.allele_specific_binding_thresholds,
        args.allele_specific_anchors,
        args.anchor_contribution_threshold,
        'pVACseq'
    ).execute()

if __name__ == "__main__":
    main()
