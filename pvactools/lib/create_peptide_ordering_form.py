import sys
import argparse
import os
import shutil

from pvactools.tools.pvacseq.generate_protein_fasta import PvacseqGenerateProteinFasta
from pvactools.lib.calculate_manufacturability import CalculateManufacturability
from pvactools.lib.generate_reviews_files import main as run_generate_reviews_files
from pvactools.lib.color_peptides51mer import main as run_color_peptides
from pvactools.lib.run_argument_utils import downstream_sequence_length, aggregate_report_evaluations

class CreatePeptideOrderingForm:
    def __init__(self, **kwargs):
        self.input_vcf = kwargs['input_vcf']
        self.flanking_sequence_length = kwargs['flanking_sequence_length']
        self.classI_aggregated_tsv = kwargs['classI_aggregated_tsv']
        self.classII_aggregated_tsv = kwargs['classII_aggregated_tsv']
        self.output_file_prefix = kwargs['output_file_prefix']
        self.sample_name = kwargs['sample_name']
        self.output_path = kwargs.pop('output_path', None)
        if self.output_path is None:
            self.output_path = f"{self.output_file_prefix}_results"
            if os.path.exists(output_path):
                if not os.path.isdir(output_path):
                    sys.exit(f"Error: {output_path} must specify a directory.")
            else:
                os.makedirs(output_path)
        self.phased_proximal_variants_vcf = kwargs.pop('phased_proximal_variants_vcf', None)
        self.external_vcf = kwargs.pop('external_vcf', None)
        self.pass_only = kwargs.pop('pass_only', False)
        self.biotypes = kwargs.pop('biotypes', ['protein_coding'])
        self.allow_incomplete_transcripts = kwargs.pop('allow_incomplete_transcripts', False)
        self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', 1000)
        self.aggregate_report_evaluation = kwargs.pop('aggregate_report_evaluation', ['Accept'])
        self.classI_IC50 = kwargs.pop('classI_IC50', 1000.0)
        self.classI_percent = kwargs.pop('classI_percent', 2.0)
        self.classII_IC50 = kwargs.pop('classII_IC50', 500.0)
        self.classII_percent = kwargs.pop('classII_percent', 2.0)
        self.prob_pos = kwargs.pop('prob_pos', [])

        self.fasta_output_file = os.path.join(self.output_path, f'{self.output_file_prefix}_{self.sample_name}.fa')
        self.combined_fasta_output_file = os.path.join(self.output_path, f'{self.output_file_prefix}_{self.sample_name}_combined.fa')
        self.peptide_manufacture_output_file = os.path.join(self.output_path, f'{self.output_file_prefix}_{self.sample_name}.manufacturability.tsv')

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            f"{tool} create_peptide_ordering_form",
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

    def execute(self):
        self.create_fastas()

        CalculateManufacturability(self.combined_fasta_output_file, self.peptide_manufacture_output_file, 'fasta').execute()

        peptide_51mer_path = run_generate_reviews_files(
            peptides_path=self.peptide_manufacture_output_file,
            classI_path=self.classI_aggregated_tsv,
            classII_path=self.classII_aggregated_tsv,
            input_vcf=self.input_vcf,
            external_vcf=self.external_vcf,
            sample_name=self.sample_name,
            allowed_evaluations=self.aggregate_report_evaluation,
            output_file_prefix=self.output_file_prefix,
            output_path=self.output_path
        )

        run_color_peptides(
            fasta_path=self.combined_fasta_output_file,
            peptides_path=peptide_51mer_path,
            sample_name=self.sample_name,
            classI_ic50_score_max=self.classI_IC50,
            classII_ic50_percentile_max=self.classI_percent,
            classII_ic50_score_max=self.classII_IC50,
            classI_ic50_percentile_max=self.classII_percent,
            problematic_position=self.prob_pos,
            output_file_prefix=self.output_file_prefix,
            output_path=self.output_path
        )

        os.remove(peptide_51mer_path)
        os.remove(self.combined_fasta_output_file)

class PvacseqCreatePeptideOrderingForm(CreatePeptideOrderingForm):
    def create_fastas(self):
        params = {
            'input_vcf': self.input_vcf,
            'sample_name': self.sample_name,
            'pass_only': self.pass_only,
            'phased_proximal_variants_vcf': self.phased_proximal_variants_vcf,
            'biotypes': self.biotypes,
            'allow_incomplete_transcripts': self.allow_incomplete_transcripts,
            'downstream_sequence_length': self.downstream_sequence_length,
            'flanking_sequence_length': self.flanking_sequence_length,
            'mutant_only': True,
            'aggregate_report_evaluation': self.aggregate_report_evaluation,
            'input_tsv': self.classI_aggregated_tsv,
            'output_file': self.fasta_output_file,
        }
        generator = PvacseqGenerateProteinFasta(**params)
        generator.generate_fasta()
        generator.trim_sequences()
        generator.filter_fasta()
        shutil.copy(generator.filtered_fasta_file_path, self.fasta_output_file)

        generator.mutant_only = False
        generator.output_file = self.combined_fasta_output_file
        generator.filter_fasta()
        shutil.copy(generator.filtered_fasta_file_path, self.combined_fasta_output_file)
        shutil.rmtree(generator.temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
