import sys

from pvactools.lib.run_argument_parser import PvacseqRunArgumentParser
from pvactools.lib.pvacseq_run_pipeline import PvacseqRunPipeline

def define_parser():
    return PvacseqRunArgumentParser().parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacseqRunPipeline(**vars(args)).execute()

if __name__ == '__main__':
    main()
