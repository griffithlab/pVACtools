import argparse
import sys

from pvactools.lib.filter import Filter, FilterCriterion

def define_parser():
    parser = argparse.ArgumentParser(
        'pvacseq transcript_support_level_filter',
        description="Filter variants processed by IEDB by transcript support level",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'input_file',
        help="The all_epitopes.tsv or filtered.tsv pVACseq report file to filter."
    )
    parser.add_argument(
        'output_file',
        help="Output .tsv file containting list of of filtered epitopes based on transcript support level."
    )
    parser.add_argument(
        "--maximum-transcript-support-level", type=int,
        help="The threshold to use for filtering epitopes on the transcript support level. "
        +"Keep all epitopes with a transcript support level <= to this cutoff.",
        default=1,
        choices=[1,2,3,4,5]
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    filter_criteria = [FilterCriterion('Transcript Support Level', '<=', args.maximum_transcript_support_level, exclude_nas=True, skip_value='Not Supported')]
    Filter(
        args.input_file,
        args.output_file,
        filter_criteria,
        ['Transcript Support Level'],
    ).execute()

if __name__ == "__main__":
    main()
