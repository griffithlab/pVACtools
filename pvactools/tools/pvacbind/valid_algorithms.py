import sys

from pvactools.lib.valid_algorithms import ValidAlgorithms

def define_parser():
    return ValidAlgorithms.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    ValidAlgorithms(args.allele, args.species).print_valid_algorithms()

if __name__ == "__main__":
    main()
