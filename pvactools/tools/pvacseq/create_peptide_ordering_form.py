import sys
import argparse
import os

from pvactools.tools.pvacseq.generate_protein_fasta import run_generate_protein_fasta
from pvactools.lib.generate_reviews_files import main as run_generate_reviews_files
from pvactools.lib.color_peptides51mer import main as run_color_peptides

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacseq create_peptide_ordering_form",
        description="Generate a peptide ordering form with coloring.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_vcf",
        help="A VEP-annotated single- or multi-sample VCF containing genotype, transcript, "
            +"Wildtype protein sequence, and Frameshift protein sequence information."
            +"The VCF may be gzipped (requires tabix index)."
    )
    parser.add_argument(
        'classI_tsv',
        help='The path to the classI all_epitopes.aggregated.tsv file with the Evaluation column filled in to mark candidates to process as Accepted'
    )
    parser.add_argument(
        'classII_tsv',
        help='The path to the classII all_epitopes.aggregated.tsv'
    )
    parser.add_argument(
        "output_file",
        help="The prefix for the output files' names"
    )
    parser.add_argument(
        "sample_name",
        help="The name of the sample being processed. Must be a sample ID in the input VCF #CHROM header line."
    )
    parser.add_argument(
        "-o", "--output-path",
        help="The path where the output will be generated. A directory will be created if not specified."
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
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence."
    )
    parser.add_argument(
        '--include-review-candidates',
        help="Include the processing of candidates marked for review.",
        default=False,
        action='store_true'
    )
    parser.add_argument(
        '--all-epitopes',
        help='If you want to generate 51mer for all epitopes',
        default=False,
        action='store_true'
    )
    parser.add_argument(
        '--classI-IC50',
        help='Maximum classI IC50 score to annotate',
        default=1000,
        type=float
    )
    parser.add_argument(
        '--classI-percent',
        help='Maximum classI percentile to annotate',
        default=2,
        type=float
    )
    parser.add_argument(
        '--classII-IC50',
        help='Maximum classII IC50 score to annotate',
        default=500,
        type=float
    )
    parser.add_argument(
        '--classII-percent',
        help='Maximum classII percentile to annotate',
        default=2,
        type=float
    )
    parser.add_argument(
        '--prob-pos',
        nargs='*', 
        help='Problematic position to make large',
        default=''
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)
    flanking_sequence_length = 25

    if not args.output_path:
        output_dir = f"{args.output_file}_results"
        os.makedirs(output_dir, exist_ok=True)
        output_path = output_dir
    else:
        output_path = args.output_path

    run_generate_protein_fasta(
        input_vcf=args.input_vcf,
        flanking_sequence_length=flanking_sequence_length,
        output_file=os.path.join(output_path, args.output_file),
        input_tsv=args.classI_tsv,
        phased_proximal_variants_vcf=args.phased_proximal_variants_vcf,
        pass_only=args.pass_only,
        biotypes=args.biotypes,
        mutant_only=True,
        aggregate_report_evaluation=['Accept', 'Review'] if args.include_review_candidates else ['Accept'],
        downstream_sequence_length=args.downstream_sequence_length,
        sample_name=args.sample_name,
        peptide_ordering_form=True
    )

    file_path = os.path.join(output_path, args.output_file)
    file_prefix = os.path.join(output_path, f"{args.output_file}_{args.sample_name}")

    os.rename(file_path, f"{file_prefix}.fa")
    os.rename(f"{file_path}_combined", f"{file_prefix}_combined.fa")
    os.rename(f"{file_path}.manufacturability.tsv", f"{file_prefix}.manufacturability.tsv")
    peptide_manufacture_path = f"{file_prefix}.manufacturability.tsv"
    combined_fasta_path = f"{file_prefix}_combined.fa"

    peptide_51mer_path = run_generate_reviews_files(
        peptides_path=peptide_manufacture_path,
        classI_path=args.classI_tsv,
        classII_path=args.classII_tsv,
        sample_name=args.sample_name,
        all_epitopes_flag=args.all_epitopes
    )

    run_color_peptides(
        fasta_path=combined_fasta_path,
        peptides_path=peptide_51mer_path,
        sample_name=args.sample_name,
        classI_ic50_score_max=args.classI_IC50,
        classII_ic50_percentile_max=args.classI_percent,
        classII_ic50_score_max=args.classII_IC50,
        classI_ic50_percentile_max=args.classII_percent,
        problematic_position=args.prob_pos,
        output_file=args.output_file,
        output_path=output_path
    )

    os.remove(peptide_51mer_path)
    os.remove(combined_fasta_path)
    

if __name__ == "__main__":
    main()