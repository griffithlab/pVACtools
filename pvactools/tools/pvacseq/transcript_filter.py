import argparse
import sys

from pvactools.lib.transcript_filter import TranscriptFilter

def define_parser():
    return TranscriptFilter.parser('pvacseq')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    TranscriptFilter(
        args.input_file,
        args.output_file,
        args.transcript_prioritization_strategy,
        args.maximum_transcript_support_level
    ).execute()

if __name__ == "__main__":
    main()
