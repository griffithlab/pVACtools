import sys
import argparse

from pvactools.lib.generate_transcripts_fasta import PvacfuseGenerateTranscriptsFasta
from pvactools.lib.run_argument_utils import downstream_sequence_length

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacfuse generate_transcripts_fasta",
        description="Generate a fasta file of matched fusion and 5'/3' wildtype transcript sequences from AGFusion or Arriba output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input",
        help="An AGFusion output directory or Arriba fusion.tsv output file."
    )
    parser.add_argument(
        "ref_fasta",
        help="A reference CDS FASTA file. Note: this input should match the build and Ensembl version used to create the fusion annotations."
    )
    parser.add_argument(
        "output_file",
        help="The output fasta file."
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshift fusion when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence.",
        type=downstream_sequence_length()
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    params = {
        'input': args.input,
        'ref_fasta': args.ref_fasta,
        'downstream_sequence_length': args.downstream_sequence_length,
        'output_file': args.output_file,
    }
    PvacfuseGenerateTranscriptsFasta(**params).execute()

if __name__ == '__main__':
    main()
