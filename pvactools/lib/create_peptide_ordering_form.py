import sys
import argparse
import os
import shutil

from pvactools.lib.generate_protein_fasta import PvacseqGenerateProteinFasta, PvacspliceGenerateProteinFasta, PvacfuseGenerateProteinFasta
from pvactools.lib.generate_reviews_files import main as run_generate_reviews_files
from pvactools.lib.color_peptides51mer import main as run_color_peptides
from pvactools.lib.run_argument_utils import downstream_sequence_length, aggregate_report_evaluations, pvacsplice_anchors

class CreatePeptideOrderingForm:
    def __init__(self, **kwargs):
        self.flanking_sequence_length = kwargs['flanking_sequence_length']
        self.classI_aggregated_tsv = kwargs['classI_aggregated_tsv']
        self.classII_aggregated_tsv = kwargs['classII_aggregated_tsv']
        self.output_file_prefix = kwargs['output_file_prefix']
        self.sample_name = kwargs['sample_name']
        self.output_path = kwargs.pop('output_path', None)
        if self.output_path is None:
            self.output_path = f"{self.output_file_prefix}_results"
            if os.path.exists(self.output_path):
                if not os.path.isdir(self.output_path):
                    sys.exit(f"Error: {self.output_path} must specify a directory.")
            else:
                os.makedirs(self.output_path)
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
        if tool == 'pvacseq':
            parser.add_argument(
                "input_vcf",
                help="A VEP-annotated single- or multi-sample VCF containing genotype, transcript, "
                    +"Wildtype protein sequence, and Frameshift protein sequence information. "
                    +"The VCF may be gzipped (requires tabix index). This VCF will be used to extract "
                    +"peptide sequences for processable variants with 25 flanking amino acids on either "
                    +"side of the mutation. These sequences will be included in the peptide ordering spreadsheet."
            )
        elif tool == 'pvacsplice':
            parser.add_argument(
                "input_file",
                help="RegTools junctions output TSV file"
            )
            parser.add_argument(
                "annotated_vcf",
                help="A VEP-annotated single- or multi-sample VCF containing genotype and transcript information."
                + "The VCF ma be gzipped (requires tabix index)."
            )
            parser.add_argument(
                "ref_fasta",
                help="A reference FASTA file. Note: this input should be the same as the RegTools vcf input."
            )
            parser.add_argument(
                "gtf_file",
                help="A reference GTF file. Note: this input should be the same as the RegTools gtf input."
            )
        elif tool == 'pvacfuse':
            parser.add_argument(
                "input",
                help="An AGFusion output directory or Arriba fusion.tsv output file."
            )
            parser.add_argument(
                "ref_fasta",
                help="A reference CDS FASTA file. Note: this input should match the build and Ensembl version used to create the fusion annotations."
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
        if tool == 'pvacseq':
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
        if tool in ['pvacseq', 'pvacsplice']:
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
        if tool == 'pvacsplice':
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
    def __init__(self, **kwargs):
        self.input_vcf = kwargs['input_vcf']
        super().__init__(**kwargs)

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
        generator.execute()
        shutil.copy(generator.manufacturability_file, self.peptide_manufacture_output_file)
        os.remove(generator.manufacturability_file)

        params['mutant_only'] = False
        params['output_file'] = self.combined_fasta_output_file
        generator = PvacseqGenerateProteinFasta(**params)
        generator.execute()
        os.remove(generator.manufacturability_file)

class PvacspliceCreatePeptideOrderingForm(CreatePeptideOrderingForm):
    def __init__(self, **kwargs):
        self.input_file = kwargs['input_file']
        self.annotated_vcf = kwargs['annotated_vcf']
        self.ref_fasta = kwargs['ref_fasta']
        self.gtf_file = kwargs['gtf_file']
        self.junction_score = kwargs.pop('junction_score', 10)
        self.variant_distance = kwargs.pop('variant_distance', 100)
        self.anchor_types = kwargs.pop('anchor_types', ['A', 'D', 'NDA'])
        self.input_vcf = None
        super().__init__(**kwargs)

    def create_fastas(self):
        params = {
            'input_file': self.input_file,
            'annotated_vcf': self.annotated_vcf,
            'ref_fasta': self.ref_fasta,
            'gtf_file': self.gtf_file,
            'junction_score': self.junction_score,
            'variant_distance': self.variant_distance,
            'anchor_types': self.anchor_types,
            'sample_name': self.sample_name,
            'pass_only': self.pass_only,
            'biotypes': self.biotypes,
            'allow_incomplete_transcripts': self.allow_incomplete_transcripts,
            'downstream_sequence_length': self.downstream_sequence_length,
            'flanking_sequence_length': self.flanking_sequence_length,
            'mutant_only': True,
            'aggregate_report_evaluation': self.aggregate_report_evaluation,
            'input_tsv': self.classI_aggregated_tsv,
            'output_file': self.fasta_output_file,
        }
        generator = PvacspliceGenerateProteinFasta(**params)
        generator.execute()
        shutil.copy(generator.manufacturability_file, self.peptide_manufacture_output_file)
        os.remove(generator.manufacturability_file)

        params['mutant_only'] = False
        params['output_file'] = self.combined_fasta_output_file
        generator = PvacspliceGenerateProteinFasta(**params)
        generator.execute()
        os.remove(generator.manufacturability_file)

class PvacfuseCreatePeptideOrderingForm(CreatePeptideOrderingForm):
    def __init__(self, **kwargs):
        self.input = kwargs['input']
        self.ref_fasta = kwargs['ref_fasta']
        self.input_vcf = None
        super().__init__(**kwargs)

    def create_fastas(self):
        params = {
            'input': self.input,
            'ref_fasta': self.ref_fasta,
            'sample_name': self.sample_name,
            'downstream_sequence_length': self.downstream_sequence_length,
            'flanking_sequence_length': self.flanking_sequence_length,
            'mutant_only': True,
            'aggregate_report_evaluation': self.aggregate_report_evaluation,
            'input_tsv': self.classI_aggregated_tsv,
            'output_file': self.fasta_output_file,
        }
        generator = PvacfuseGenerateProteinFasta(**params)
        generator.execute()
        shutil.copy(generator.manufacturability_file, self.peptide_manufacture_output_file)
        os.remove(generator.manufacturability_file)

        params['mutant_only'] = False
        params['output_file'] = self.combined_fasta_output_file
        generator = PvacfuseGenerateProteinFasta(**params)
        generator.execute()
        os.remove(generator.manufacturability_file)

if __name__ == "__main__":
    main()
