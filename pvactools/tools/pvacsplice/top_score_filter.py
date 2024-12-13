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
        maximum_transcript_support_level=args.maximum_transcript_support_level,
    ).execute()

if __name__ == "__main__":
    main()
