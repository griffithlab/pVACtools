import sys
import argparse

from pvactools.lib.generate_protein_fasta import PvacfuseGenerateProteinFasta
from pvactools.lib.run_argument_utils import aggregate_report_evaluations

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacfuse generate_protein_fasta",
        description="Generate a fasta file with a specific flanking sequence length around the fusion",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "transcripts_fasta",
        help="A pVACfuse transcripts.fa file with transcript protein sequences of fusions and matching 5'/3' wildtypes. "
             + "This file can be found in the top-level output directory of your pVACfuse run can be generated using the `pvacseq generate_transcripts_fasta` command."
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
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    params = {
        'transcripts_fasta': args.transcripts_fasta,
        'flanking_sequence_length': args.flanking_sequence_length,
        'mutant_only': args.mutant_only,
        'aggregate_report_evaluation': args.aggregate_report_evaluation,
        'output_file': args.output_file,
    }
    PvacfuseGenerateProteinFasta(**params).execute()

if __name__ == '__main__':
    main()
