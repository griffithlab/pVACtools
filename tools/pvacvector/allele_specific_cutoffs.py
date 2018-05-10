from lib.allele_specific_cutoffs import *

def define_parser():
    return AlleleSpecificCutoffs.parser('pvacvector')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    AlleleSpecificCutoffs(args.allele).print_allele_specific_cutoffs()

if __name__ == "__main__":
    main()
