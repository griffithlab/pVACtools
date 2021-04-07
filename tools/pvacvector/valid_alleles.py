from lib.valid_alleles import *

def define_parser():
    return ValidAlleles.parser('pvacvector')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    ValidAlleles(args.prediction_algorithm, args.species).print_valid_alleles()

if __name__ == "__main__":
    main()
