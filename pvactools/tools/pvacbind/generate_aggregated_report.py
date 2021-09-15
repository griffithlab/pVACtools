import sys
import argparse
import tempfile

from pvactools.lib.aggregate_all_epitopes import UnmatchedSequenceAggregateAllEpitopes

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacbind generate_aggregated_report",
        description="Generate an aggregated report from a pVACbind .all_epitopes.tsv report file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input_file",
        help="A pVACbind .all_epitopes.tsv report file"
    )
    parser.add_argument(
        "output_file",
        help="The file path to write the aggregated report tsv to"
    )

    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    tmp_fh = tempfile.NamedTemporaryFile()

    print("Creating Aggreggated Report")
    UnmatchedSequenceAggregateAllEpitopes(args.input_file, args.output_file).execute()
    print("Completed")

if __name__ == '__main__':
    main()
