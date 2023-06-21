import sys
import argparse
import tempfile

from pvactools.lib.identify_problematic_amino_acids import IdentifyProblematicAminoAcids

def define_parser():
    return IdentifyProblematicAminoAcids.parser('pvacbind')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    IdentifyProblematicAminoAcids(args.input_file, args.output_file, args.problematic_amino_acids, 'pVACbind').execute()

if __name__ == "__main__":
    main()
