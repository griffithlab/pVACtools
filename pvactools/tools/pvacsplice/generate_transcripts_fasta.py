import sys
import argparse

from pvactools.lib.generate_transcripts_fasta import PvacspliceGenerateTranscriptsFasta
from pvactools.lib.run_argument_utils import pvacsplice_anchors, downstream_sequence_length

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacsplice generate_transcripts_fasta",
        description="Generate afasta file from a RegTools junctions output TSV file with transcripts sequences of splicing events and matching wildtypes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_file",
        help="RegTools junctions output TSV file"
    )
    parser.add_argument(
        "output_file",
        help="The output fasta file."
    )
    parser.add_argument(
        "annotated_vcf",
        help="A VEP-annotated single- or multi-sample VCF containing genotype and transcript information."
        + "The VCF may be gzipped (requires tabix index)."
    )
    parser.add_argument(
        "ref_fasta",
        help="A reference FASTA file. Note: this input should be the same as the RegTools vcf input."
    )
    parser.add_argument(
        "gtf_file",
        help="A reference GTF file. Note: this input should be the same as the RegTools gtf input."
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
        "-j", "--junction-score", type=int,
        help="Junction Coverage Cutoff. Only sites above this read depth cutoff will be considered.",
        default=10
    )
    parser.add_argument(
        "-v", "--variant-distance", type=int,
        help="Regulatory variants can lie inside or outside of splicing junction."
        + "Maximum distance window (upstream and downstream) for a variant outside the junction.",
        default=100
    )
    parser.add_argument(
        "--anchor-types", type=pvacsplice_anchors(),
        help="The anchor types of junctions to use. Multiple anchors can be specified using a comma-separated list."
        + "Choices: A, D, NDA, DA, N",
        default=['A', 'D', 'NDA'],
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshift splice sites when creating the fasta file. "
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
        'input_file': args.input_file,
        'annotated_vcf': args.annotated_vcf,
        'ref_fasta': args.ref_fasta,
        'gtf_file': args.gtf_file,
        'sample_name': args.sample_name,
        'pass_only': args.pass_only,
        'biotypes': args.biotypes,
        'allow_incomplete_transcripts': args.allow_incomplete_transcripts,
        'downstream_sequence_length': args.downstream_sequence_length,
        'junction_score': args.junction_score,
        'variant_distance': args.variant_distance,
        'anchor_types': args.anchor_types,
        'output_file': args.output_file,
    }
    PvacspliceGenerateTranscriptsFasta(**params).execute()

if __name__ == '__main__':
    main()
