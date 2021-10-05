import sys

from pvactools.lib.allele_specific_cutoffs import AlleleSpecificCutoffs

def define_parser():
    return AlleleSpecificCutoffs.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    AlleleSpecificCutoffs(args.allele).print_allele_specific_cutoffs()

if __name__ == "__main__":
    main()
