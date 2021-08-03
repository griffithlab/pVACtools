import sys
from lib.net_chop import *

def define_parser():
    return NetChop.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    NetChop(args.input_file, args.output_file, args.method, args.threshold).execute()

if __name__ == "__main__":
    main()