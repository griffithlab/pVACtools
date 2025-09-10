import sys
import argparse
import os

from pvactools.tools.pvacseq.generate_protein_fasta import run_generate_protein_fasta
from pvactools.lib.generate_reviews_files import main as run_generate_reviews_files
from pvactools.lib.color_peptides51mer import main as run_color_peptides
from pvactools.lib.run_utils import aggregate_report_evaluations

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacseq create_peptide_ordering_form",
        description="Generate peptide ordering files (FASTA, annotated ordering Excel spreadsheet, and review template Excel spreadsheet) to streamline preparation of peptides for synthesis and review.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_vcf",
        help="A VEP-annotated single- or multi-sample VCF containing genotype, transcript, "
            +"Wildtype protein sequence, and Frameshift protein sequence information. "
            +"The VCF may be gzipped (requires tabix index). This VCF will be used to extract "
            +"peptide sequences for processable variants with 25 flanking amino acids on either "
            +"side of the mutation. These sequences will be included in the peptide ordering spreadsheet."
    )
    parser.add_argument(
        "flanking_sequence_length",
        help="Number of amino acids to add on each side of the mutation when creating the FASTA.",
        type=int,
    )
    parser.add_argument(
        'classI_aggregated_tsv',
        help="The path to the classI all_epitopes.aggregated.tsv file with the Evaluation column filled in to mark candidates "
            +"to process as 'Accept'. Only candidates marked as Accept in this file will be included in the ordering "
            +"spreadsheet. This file is commonly created by importing the aggregated class I report from pVACseq into pVACview, "
            +"investigating candidates, selecting appropriate evaluations, and exporting the results in TSV format."
    )
    parser.add_argument(
        'classII_aggregated_tsv',
        help='The path to the classII all_epitopes.aggregated.tsv'
    )
    parser.add_argument(
        "output_file_prefix",
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
        help="A VCF with phased proximal variant information to incorporate into the predicted fasta sequences "
            +"generated from the input_vcf. Must be gzipped and tabix indexed."
    )
    parser.add_argument(
        '--external-vcf',
        help='A VCF file from an external provider to check variants against. Any variant '
            +'with a PASS filter or no other filter applied will be marked as called in the '
            +'"Variant Called in External VCF" column of the updated aggregated report '
            +'"<sample_name>.Annotated.Neoantigen_Candidates.xlsx"'
    )
    parser.add_argument(
        '--pass-only',
        help="Only process VCF entries with a PASS status.",
        default=False,
        action='store_true',
    )
    parser.add_argument(
        "--biotypes", type=lambda s:[a for a in s.split(',')],
        help="A list of biotypes to use for pre-filtering transcripts when generating peptide sequences from "
            +"the input_vcf.",
        default=['protein_coding']
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence."
    )
    parser.add_argument(
        "--aggregate-report-evaluation",
        help="Only include variants where the Evaluation column in the classI_aggregated_tsv matches this evaluation. "
            +"Valid values for this field are Accept, Reject, Pending, and Review. Specify multiple values as "
            +"a comma-separated list to include multiple evaluation states.",
        default='Accept',
        type=aggregate_report_evaluations(),
    )
    parser.add_argument(
        '--classI-IC50',
        help="Bold the Best Peptide from the classI_aggregated_tsv file in the 'CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES' "
            +"column of the ordering spreadsheet only if the IC50 score is less than this cutoff or the --classI-percent cutoff is met.",
        default=1000,
        type=float
    )
    parser.add_argument(
        '--classI-percent',
        help="Color the Best Peptide from the classI_aggregated_tsv file in the 'CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES' "
            +"column of the ordering spreadsheet only if this percentile cutoff is met or the IC50 score is below the specified --classI-IC50 maximum.",
        default=2,
        type=float
    )
    parser.add_argument(
        '--classII-IC50',
        help="Bold the Best Peptide from the classII_aggregated_tsv file in the 'CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES' "
            +"column of the ordering spreadsheet only if the IC50 score is less than this cutoff or the --classII-percent cutoff is met.",
        default=500,
        type=float
    )
    parser.add_argument(
        '--classII-percent',
        help="Bold the Best Peptide from the classII_aggregated_tsv file in the 'CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES' "
            +"column of the ordering spreadsheet only if this percentile cutoff is met or the IC50 score is below the specified --classII-IC50 maximum.",
        default=2,
        type=float
    )
    parser.add_argument(
        '--prob-pos',
        type=lambda s: [item.strip() for item in s.split(',')],
        help='Comma-separated list of problematic positions to make large in the ordering spreadsheet.',
        default=[]
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    output_path = args.output_path if args.output_path else f"{args.output_file_prefix}_results"
    if os.path.exists(output_path):
        if not os.path.isdir(output_path):
            sys.exit(f"Error: {output_path} must specify a directory.")
    else:
        os.makedirs(output_path)

    run_generate_protein_fasta(
        input_vcf=args.input_vcf,
        flanking_sequence_length=args.flanking_sequence_length,
        output_file=os.path.join(output_path, args.output_file_prefix),
        input_tsv=args.classI_aggregated_tsv,
        phased_proximal_variants_vcf=args.phased_proximal_variants_vcf,
        pass_only=args.pass_only,
        biotypes=args.biotypes,
        mutant_only=True,
        aggregate_report_evaluation=args.aggregate_report_evaluation,
        downstream_sequence_length=args.downstream_sequence_length,
        sample_name=args.sample_name,
        peptide_ordering_form=True
    )

    file_path = os.path.join(output_path, args.output_file_prefix)
    file_prefix = os.path.join(output_path, f"{args.output_file_prefix}_{args.sample_name}")

    os.rename(file_path, f"{file_prefix}.fa")
    os.rename(f"{file_path}_combined", f"{file_prefix}_combined.fa")
    os.rename(f"{file_path}.manufacturability.tsv", f"{file_prefix}.manufacturability.tsv")
    peptide_manufacture_path = f"{file_prefix}.manufacturability.tsv"
    combined_fasta_path = f"{file_prefix}_combined.fa"

    peptide_51mer_path = run_generate_reviews_files(
        peptides_path=peptide_manufacture_path,
        classI_path=args.classI_aggregated_tsv,
        classII_path=args.classII_aggregated_tsv,
        input_vcf=args.input_vcf,
        external_vcf=args.external_vcf,
        sample_name=args.sample_name,
        allowed_evaluations=args.aggregate_report_evaluation,
        output_file_prefix=args.output_file_prefix,
        output_path=output_path
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
        output_file_prefix=args.output_file_prefix,
        output_path=output_path
    )

    os.remove(peptide_51mer_path)
    os.remove(combined_fasta_path)
    

if __name__ == "__main__":
    main()