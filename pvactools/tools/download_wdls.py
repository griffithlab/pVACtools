import argparse
import wget
import os
import sys

def define_parser():
    parser = argparse.ArgumentParser(
        "pvactools download_wdls",
        description="Download pVACtools WDLs to run the main pVACseq and pVACfuse pipelines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('destination_directory', help='Directory for downloading WDLs')
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    final_dir = os.path.join(args.destination_directory, 'pvactools_wdls')
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)

    wget.download('https://raw.githubusercontent.com/wustl-oncology/analysis-wdls/master/definitions/tools/pvacseq.wdl', final_dir)
    wget.download('https://raw.githubusercontent.com/wustl-oncology/analysis-wdls/master/definitions/tools/pvacfuse.wdl', final_dir)
    print("")

if __name__ == '__main__':
    main()
