import sys

from pvactools.lib.top_score_filter import PvacbindTopScoreFilter, TopScoreFilter

def define_parser():
    return TopScoreFilter.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacbindTopScoreFilter(args.input_file, args.output_file, top_score_metric=args.top_score_metric, top_score_metric2=args.top_score_metric2).execute()

if __name__ == "__main__":
    main()
