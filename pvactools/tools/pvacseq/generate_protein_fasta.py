import sys
import argparse

from pvactools.lib.generate_protein_fasta import PvacseqGenerateProteinFasta
from pvactools.lib.run_utils import *

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacseq generate_protein_fasta",
        description="Generate an annotated fasta file from a VCF with protein sequences of mutations and matching wildtypes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input_vcf",
        help="A VEP-annotated single- or multi-sample VCF containing genotype, transcript, "
            +"Wildtype protein sequence, and Frameshift protein sequence information."
            +"The VCF may be gzipped (requires tabix index)."
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
        help = "A pVACseq all_epitopes, filtered, or aggregated TSV file with epitopes to use for subsetting the input VCF to peptides of interest. Only the peptide sequences for the epitopes in the TSV will be used when creating the FASTA. When running with an aggregated TSV, the sequences will be further narrowed down to only include variants with the selected --aggregate-report-evaluation."
    )
    parser.add_argument(
        "-p", "--phased-proximal-variants-vcf",
        help="A VCF with phased proximal variant information to incorporate into the predicted fasta sequences. Must be gzipped and tabix indexed."
    )
    parser.add_argument(
        '--pass-only',
        help="Only process VCF entries with a PASS status.",
        default=False,
        action='store_true',
    )
    parser.add_argument(
        "--biotypes", type=lambda s:[a for a in s.split(',')],
        help="A list of biotypes to use for pre-filtering transcripts for processing in the pipeline.",
        default=['protein_coding']
    )
    parser.add_argument(
        "--allow-incomplete-transcripts",
        help="By default, transcripts annotated with incomplete CDS (i.e., 'cds_start_NF' or 'cds_end_NF' flags in the VEP CSQ field) "
                + "are excluded from analysis, as they often produce invalid protein sequences. "
                + "Use this flag to allow candidates from such transcripts. Only peptides that do not contain 'X' will be included. "
                + "These candidates will be deprioritized relative to those from transcripts without incomplete CDS flags.",
        default=False,
        action='store_true'
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
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence."
    )
    parser.add_argument(
        "-s", "--sample-name",
        help="The name of the sample being processed. Required when processing a multi-sample VCF and must be a sample ID in the input VCF #CHROM header line."
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
        'input_vcf': args.input_vcf,
        'sample_name': args.sample_name,
        'pass_only': args.pass_only,
        'phased_proximal_variants_vcf': args.phased_proximal_variants_vcf,
        'biotypes': args.biotypes,
        'allow_incomplete_transcripts': args.allow_incomplete_transcripts,
        'downstream_sequence_length': downstream_sequence_length,
        'flanking_sequence_length': args.flanking_sequence_length,
        'mutant_only': args.mutant_only,
        'aggregate_report_evaluation': args.aggregate_report_evaluation,
        'input_tsv': args.input_tsv,
        'output_file': args.output_file,
    }
    PvacseqGenerateProteinFasta(**params).execute()


if __name__ == '__main__':
    main()
