import argparse
import wget
import os
import sys

def define_parser():
    parser = argparse.ArgumentParser("pvactools download_cwls")
    parser.add_argument('destination_directory', help='Directory for downloading CWLs')
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    final_dir = os.path.join(args.destination_directory, 'pvactools_cwls')
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)

    wget.download('https://raw.githubusercontent.com/genome/cancer-genomics-workflow/master/pvactools/pvacseq.cwl', final_dir)
    wget.download('https://raw.githubusercontent.com/genome/cancer-genomics-workflow/master/pvactools/pvacfuse.cwl', final_dir)
    wget.download('https://raw.githubusercontent.com/genome/cancer-genomics-workflow/master/pvactools/pvacvector.cwl', final_dir)
    print("")

if __name__ == '__main__':
    main()
