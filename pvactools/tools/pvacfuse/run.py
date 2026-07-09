import sys

from pvactools.lib.run_argument_parser import PvacfuseRunArgumentParser
from pvactools.lib.pvacfuse_run_pipeline import PvacfuseRunPipeline

def define_parser():
    return PvacfuseRunArgumentParser().parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    PvacfuseRunPipeline(**vars(args)).execute()

if __name__ == '__main__':
    main()
