import sys
import argparse

from pvactools.lib.generate_transcripts_fasta import PvacseqGenerateTranscriptsFasta
from pvactools.lib.run_argument_utils import downstream_sequence_length

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacseq generate_transcripts_fasta",
        description="Generate a fasta file from a VCF with transcript sequences of mutations and matching wildtypes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input_vcf",
        help="A VEP-annotated single- or multi-sample VCF containing genotype, transcript, "
            +"Wildtype protein sequence, and Frameshift protein sequence information."
            +"The VCF may be gzipped (requires tabix index)."
    )
    parser.add_argument(
        "output_file",
        help="The output fasta file."
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
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence.",
        type=downstream_sequence_length()
    )
    parser.add_argument(
        "-s", "--sample-name",
        help="The name of the sample being processed. Required when processing a multi-sample VCF and must be a sample ID in the input VCF #CHROM header line."
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    params = {
        'input_vcf': args.input_vcf,
        'sample_name': args.sample_name,
        'pass_only': args.pass_only,
        'phased_proximal_variants_vcf': args.phased_proximal_variants_vcf,
        'biotypes': args.biotypes,
        'allow_incomplete_transcripts': args.allow_incomplete_transcripts,
        'downstream_sequence_length': args.downstream_sequence_length,
        'output_file': args.output_file,
    }
    PvacseqGenerateTranscriptsFasta(**params).execute()


if __name__ == '__main__':
    main()
