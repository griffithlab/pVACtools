import argparse
import wget
import os
import sys

def define_parser():
    parser = argparse.ArgumentParser(
        "pvactools download_cwls",
        description="Download pVACtools CWLs for each tool's main pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('destination_directory', help='Directory for downloading CWLs')
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    final_dir = os.path.join(args.destination_directory, 'pvactools_cwls')
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)

    wget.download('https://raw.githubusercontent.com/genome/analysis-workflows/master/definitions/tools/pvacseq.cwl', final_dir)
    wget.download('https://raw.githubusercontent.com/genome/analysis-workflows/master/definitions/tools/pvacfuse.cwl', final_dir)
    wget.download('https://raw.githubusercontent.com/genome/analysis-workflows/master/definitions/tools/pvacvector.cwl', final_dir)
    print("")

if __name__ == '__main__':
    main()
