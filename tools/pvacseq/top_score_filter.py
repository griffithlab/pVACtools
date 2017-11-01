import sys
from lib.top_score_filter import *

def define_parser():
    return TopScoreFilter.parser('pvacseq')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    TopScoreFilter(args.input_file, args.output_file, args.top_score_metric).execute()

if __name__ == "__main__":
    main()
