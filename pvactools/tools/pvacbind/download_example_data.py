import sys

from pvactools.lib.download_example_data import DownloadExampleData

def define_parser():
    return DownloadExampleData.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    DownloadExampleData(args.destination_directory , 'pvacbind').execute()


if __name__ == '__main__':
    main()
