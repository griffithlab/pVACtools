import argparse
import sys
from lib.filter import *

def define_parser():
    parser = argparse.ArgumentParser('pvacseq transcript_support_level_filter')
    parser.add_argument(
        'input_file',
        help="The all_epitopes.tsv or filtered.tsv pVACseq report file to filter"
    )
    parser.add_argument(
        'output_file',
        help="Output .tsv file containting list of of filtered epitopes based on transcript support level."
    )
    parser.add_argument(
        "--maximum-transcript-support-level", type=int,
        help="The threshold to use for filtering epitopes on the transcript support level. "
        +"Keep all epitopes with a transcript support level <= to this cutoff. Default: 1",
        default=1,
        choices=[1,2,3,4,5]
    )
    parser.add_argument(
        '--exclude-NAs',
        help="Exclude NA values from the filtered output. Default: False",
        default=False,
        action='store_true'
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    filter_criteria = [{'column': 'Transcript Support Level', 'operator': '<=', 'threshold': args.maximum_transcript_support_level}]
    Filter(
        args.input_file,
        args.output_file,
        filter_criteria,
        args.exclude_NAs,
        ['Transcript Support Level'],
    ).execute()

if __name__ == "__main__":
    main()
