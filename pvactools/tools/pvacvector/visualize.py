import argparse
import os
import sys

from pvactools.lib.vector_visualization import VectorVisualization

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacvector visualize",
        description=" Create a visualization of the vector from a pVACvector result fasta",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_fasta",
        help="A pVACvector result FASTA file to visualize",
    )
    parser.add_argument(
        "output_directory",
        help="The output directory to save the visualization graphic to",
    )
    parser.add_argument(
        '-s', '--spacers', type=lambda s:[spacer for spacer in s.split(',')],
        help="Comma-separated list of peptides that are used as spacers in the pVACvector result fasta file",
        default="AAY,HHHH,GGS,GPGPG,HHAA,AAL,HH,HHC,HHH,HHHD,HHL,HHHC"
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if 'DISPLAY' in os.environ.keys():
        VectorVisualization(args.input_fasta, args.output_directory, args.spacers).draw()
    else:
        raise Exception("DISPLAY environment variable not set. Unable to create vector visualization.")

if __name__ == "__main__":
    main()
