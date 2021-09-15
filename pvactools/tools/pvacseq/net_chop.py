import sys

from pvactools.lib.net_chop import NetChop

def define_parser():
    return NetChop.parser('pvacseq')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    NetChop(args.input_file, args.input_fasta, args.output_file, args.method, args.threshold, 'pVACseq').execute()

if __name__ == "__main__":
    main()
