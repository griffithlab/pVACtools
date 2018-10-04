import sys
import argparse
import tempfile
from lib.condense_final_report import *
from lib.rank_epitopes import *

def define_parser():
    parser = argparse.ArgumentParser("pvacseq generate_condensed_ranked_report")

    parser.add_argument(
        "input_file",
        help="A pVACseq .all_epitopes.tsv or .filtered.tsv report file"
    )
    parser.add_argument(
        "output_file",
        help="The file path to write the condensed, ranked report tsv to"
    )
    parser.add_argument(
        '-m', '--top-score-metric',
        choices=['lowest', 'median'],
        default='median',
        help="The ic50 scoring metric to use for ranking epitopes by binding-threshold and minimum fold change. "
             + "lowest: Best MT Score/Corresponding Fold Change - lowest MT ic50 binding score/corresponding fold change of all chosen prediction methods. "
             + "median: Median MT Score/Median Fold Change - median MT ic50 binding score/fold change of all chosen prediction methods. "
             + "Default: median"
    )

    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    tmp_fh = tempfile.NamedTemporaryFile()

    print("Creating Condensed Report")
    CondenseFinalReport(args.input_file, tmp_fh.name, args.top_score_metric).execute()
    print("Completed")

    print("Ranking neoepitopes")
    RankEpitopes(tmp_fh.name, args.output_file).execute()
    print("Completed")

    tmp_fh.close()

if __name__ == '__main__':
    main()
