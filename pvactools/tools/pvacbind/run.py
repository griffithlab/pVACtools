import sys

from pvactools.lib.run_argument_parser import PvacbindRunArgumentParser
from pvactools.lib.pvacbind_run_pipeline import PvacbindRunPipeline

def define_parser():
    return PvacbindRunArgumentParser().parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacbindRunPipeline(**vars(args)).execute()

if __name__ == '__main__':
    main()
