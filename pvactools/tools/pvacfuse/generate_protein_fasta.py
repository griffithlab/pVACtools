import sys
import argparse

from pvactools.lib.generate_protein_fasta import PvacfuseGenerateProteinFasta
from pvactools.lib.run_utils import *

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacfuse generate_protein_fasta",
        description="Generate an annotated fasta file from AGFusion or Arriba output.",
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
        "flanking_sequence_length", type=int,
        help="Number of amino acids to add on each side of the mutation when creating the FASTA.",
    )
    parser.add_argument(
        "output_file",
        help="The output fasta file."
    )
    parser.add_argument(
        "--input-tsv",
        help = "A pVACfuse all_epitopes, filtered, or aggregated TSV file with epitopes to use for subsetting the input file to peptides of interest. Only the peptide sequences for the variants in the TSV will be used when creating the FASTA. When running with an aggregated TSV, the sequences will be further narrowed down to only include variants with the selected --aggregate-report-evaluation."
    )
    parser.add_argument(
        "--mutant-only",
        help="Only output mutant peptide sequences",
        default=False,
        action='store_true',
    )
    parser.add_argument(
        "--aggregate-report-evaluation",
        help="When running with an aggregate report input TSV, only include variants with this evaluation. Valid values for this field are Accept, Reject, Pending, and Review. Specifiy multiple values as a comma-separated list to include multiple evaluation states.",
        default='Accept',
        type=aggregate_report_evaluations(),
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshift fusion when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence."
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.downstream_sequence_length == 'full':
        downstream_sequence_length = None
    elif args.downstream_sequence_length.isdigit():
        downstream_sequence_length = int(args.downstream_sequence_length)
    else:
        sys.exit("The downstream sequence length needs to be a positive integer or 'full'")

    params = {
        'input': args.input,
        'ref_fasta': args.ref_fasta,
        'downstream_sequence_length': downstream_sequence_length,
        'flanking_sequence_length': args.flanking_sequence_length,
        'mutant_only': args.mutant_only,
        'aggregate_report_evaluation': args.aggregate_report_evaluation,
        'input_tsv': args.input_tsv,
        'output_file': args.output_file,
    }
    PvacfuseGenerateProteinFasta(**params).execute()

if __name__ == '__main__':
    main()
