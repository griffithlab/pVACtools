import sys

from pvactools.lib.mark_genes_of_interest import MarkGenesOfInterest

def define_parser():
    return MarkGenesOfInterest.parser('pvacseq')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    MarkGenesOfInterest(args.input_file, args.output_file, args.genes_of_interest_file, file_type='pVACseq').execute()

if __name__ == "__main__":
    main()
