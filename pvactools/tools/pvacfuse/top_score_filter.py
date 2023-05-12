import sys

from pvactools.lib.top_score_filter import PvacfuseTopScoreFilter, TopScoreFilter

def define_parser():
    return TopScoreFilter.parser('pvacfuse')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacfuseTopScoreFilter(args.input_file, args.output_file, top_score_metric=args.top_score_metric).execute()

if __name__ == "__main__":
    main()
