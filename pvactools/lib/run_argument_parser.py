from abc import ABCMeta
import argparse
import textwrap

from pvactools.lib.prediction_class import PredictionClass
import pvactools.lib.net_chop
from pvactools.lib.run_utils import *

class RunArgumentParser(metaclass=ABCMeta):
    def __init__(self, tool_name, input_file_help):
        parser = argparse.ArgumentParser(
            "%s run" % tool_name,
            description="Run the {} pipeline".format(tool_name.replace('vac', 'VAC')),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            "input_file",
            help=input_file_help
        )
        if tool_name == 'pvacseq':
            sample_name_help = "The name of the tumor sample being processed. When processing a multi-sample VCF the sample name must be a sample ID in the input VCF #CHROM header line."
        else:
            sample_name_help = "The name of the sample being processed. This will be used as a prefix for output files."
        parser.add_argument(
            "sample_name",
            help=sample_name_help
        )
        parser.add_argument(
            "allele", type=lambda s:[a for a in s.split(',')],
            help="Name of the allele to use for epitope prediction. "
                 + "Multiple alleles can be specified using a comma-separated list. "
                 + "For a list of available alleles, use: `{} valid_alleles`.".format(tool_name),
        )
        parser.add_argument(
            "prediction_algorithms",
            choices=PredictionClass.prediction_methods_with_all(),
            nargs="+",
            help="The epitope prediction algorithms to use. Multiple prediction algorithms can be specified, separated by spaces.",
        )
        parser.add_argument(
            "output_dir",
            help="The directory for writing all result files."
        )
        parser.add_argument(
            "--iedb-install-directory",
            help="Directory that contains the local installation of IEDB MHC I and/or MHC II."
        )
        parser.add_argument(
            "-r", "--iedb-retries",type=int,
            default=5,
            help="Number of retries when making requests to the IEDB RESTful web interface. Must be less than or equal to 100.",
        )
        parser.add_argument(
            "-k", "--keep-tmp-files",
            action='store_true',
            help="Keep intermediate output files. This might be useful for debugging purposes.",
        )
        parser.add_argument(
            "-t", "--n-threads",type=int,
            default=1,
            help="Number of threads to use for parallelizing peptide-MHC binding prediction calls.",
        )
        parser.add_argument(
            "--netmhciipan-version",
            choices=["4.3", "4.2", "4.1", "4.0"],
            default="4.1",
            help="Specify the version of NetMHCIIpan or NetMHCIIpanEL to be used during the run.",
        )
        self.parser = parser

    def epitope_args(self):
        self.parser.add_argument(
            "-e1", "--class-i-epitope-length", type=lambda s:[int(epl) for epl in s.split(',')],
            default=[8,9,10,11],
            help="Length of MHC Class I subpeptides (neoepitopes) to predict. "
                 + "Multiple epitope lengths can be specified using a comma-separated list. "
                 + "Typical epitope lengths vary between 8-15. "
                 + "Required for Class I prediction algorithms.",
        )
        self.parser.add_argument(
                "-e2", "--class-ii-epitope-length", type=lambda s:[int(epl) for epl in s.split(',')],
            default=[12,13,14,15,16,17,18],
            help="Length of MHC Class II subpeptides (neoepitopes) to predict. "
                 + "Multiple epitope lengths can be specified using a comma-separated list. "
                 + "Typical epitope lengths vary between 11-30. "
                 + "Required for Class II prediction algorithms.",
        )

    def binding_args(self, tool_name):
        self.parser.add_argument(
            "-b","--binding-threshold", type=int,
            default=500,
            help="Report only epitopes where the mutant allele has ic50 binding scores below this value.",
        )
        self.parser.add_argument(
            '--percentile-threshold', type=float_range(0.0,100.0),
            help="Report only epitopes where the mutant allele "
                 +"has a percentile rank below this value."
        )
        self.parser.add_argument(
            '--percentile-threshold-strategy',
            choices=['conservative', 'exploratory'],
            help="Specify the candidate inclusion strategy. The 'conservative' option requires a candidate to pass BOTH the binding threshold and percentile threshold (default)."
                 + " The 'exploratory' option requires a candidate to pass EITHER the binding threshold or the percentile threshold.",
            default="conservative",
        )
        self.parser.add_argument(
            '--allele-specific-binding-thresholds',
            help="Use allele-specific binding thresholds. To print the allele-specific binding thresholds run `%s allele_specific_cutoffs`. " % tool_name
                 + "If an allele does not have a special threshold value, the `--binding-threshold` value will be used.",
            default=False,
            action='store_true',
        )
        self.parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            default='median',
            help="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. "
                 + "lowest: Use the best MT Score and Corresponding Fold Change (i.e. the lowest MT ic50 binding score and corresponding fold change of all chosen prediction methods). "
                 + "median: Use the median MT Score and Median Fold Change (i.e. the  median MT ic50 binding score and fold change of all chosen prediction methods)."
        )

    def prediction_args(self):
        self.parser.add_argument(
            '--net-chop-method',
            choices=pvactools.lib.net_chop.methods,
            default=None,
            help="NetChop prediction method to use (\"cterm\" for C term 3.0, \"20s\" for 20S 3.0). C-term 3.0 is trained with publicly available MHC class I ligands and the authors believe that is performs best in predicting the boundaries of CTL epitopes. 20S is trained with in vitro degradation data.",
        )
        self.parser.add_argument(
            '--netmhc-stab',
            action='store_true',
            help="Run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes."
        )
        self.parser.add_argument(
            '--net-chop-threshold', type=float,
            default=0.5,
            help="NetChop prediction threshold (increasing the threshold results in better specificity, but worse sensitivity).",
        )
        self.parser.add_argument(
            '--problematic-amino-acids', type=lambda s:[a for a in s.split(',')],
            help=textwrap.dedent('''\
            A list of amino acids to consider as problematic. Each entry can be specified in the following format:
            `amino_acid(s)`: One or more one-letter amino acid codes. Any occurrence of this amino acid string,
                             regardless of the position in the epitope, is problematic. When specifying more than
                             one amino acid, they will need to occur together in the specified order.
            `amino_acid:position`: A one letter amino acid code, followed by a colon separator, followed by a positive
                                   integer position (one-based). The occurrence of this amino acid at the position
                                   specified is problematic., E.g. G:2 would check for a Glycine at the second position
                                   of the epitope. The N-terminus is defined as position 1.
            `amino_acid:-position`: A one letter amino acid code, followed by a colon separator, followed by a negative
                                    integer position. The occurrence of this amino acid at the specified position from
                                    the end of the epitope is problematic. E.g., G:-3 would check for a Glycine at the
                                    third position from the end of the epitope. The C-terminus is defined as position -1.''')
        )
        self.parser.add_argument(
            '--run-reference-proteome-similarity',
            action='store_true',
            help="Blast peptides against the reference proteome."
        )
        self.parser.add_argument(
            '--blastp-path',
            help="Blastp installation path.",
        )
        self.parser.add_argument(
            '--blastp-db',
            choices=['refseq_select_prot', 'refseq_protein'],
            default='refseq_select_prot',
            help="The blastp database to use.",
        )
        self.parser.add_argument(
            '--peptide-fasta',
            help="When running the reference proteome similarity step, use this reference peptide FASTA file to find matches instead of blastp."
        )
        self.parser.add_argument(
            '-a', '--additional-report-columns',
            choices=['sample_name'],
            help="Additional columns to output in the final report. If sample_name is chosen, this will add a column with the sample name in every row of the output. This can be useful if you later want to concatenate results from multiple individuals into a single file."
        )
        self.parser.add_argument(
            "-s", "--fasta-size",type=int,
            default=200,
            help="Number of FASTA entries per IEDB request. "
                 + "For some resource-intensive prediction algorithms like Pickpocket and NetMHCpan it might be helpful to reduce this number. "
                 + "Needs to be an even number.",
        )
        self.parser.add_argument(
            '--exclude-NAs',
            help="Exclude NA values from the filtered output.",
            default=False,
            action='store_true'
        )

    def pass_only_args(self):
        self.parser.add_argument(
            '--pass-only',
            help="Only process VCF entries with a PASS status.",
            default=False,
            action='store_true'
        )

    def fasta_generation(self):
        self.parser.add_argument(
            "-d", "--downstream-sequence-length",
            default='1000',
            help="Cap to limit the downstream sequence length for frameshifts when creating the FASTA file. "
                 + "Use 'full' to include the full downstream sequence."
        )

    def expression_coverage_args(self):
        self.parser.add_argument(
            '--normal-sample-name',
            help="In a multi-sample VCF, the name of the matched normal sample."
        )
        self.parser.add_argument(
            '--normal-cov', type=int,
            help="Normal Coverage Cutoff. Only sites above this read depth cutoff will be considered.",
            default=5
        )
        self.parser.add_argument(
            '--tdna-cov', type=int,
            help="Tumor DNA Coverage Cutoff. Only sites above this read depth cutoff will be considered.",
            default=10
        )
        self.parser.add_argument(
            '--trna-cov', type=int,
            help="Tumor RNA Coverage Cutoff. Only sites above this read depth cutoff will be considered.",
            default=10
        )
        self.parser.add_argument(
            '--normal-vaf', type=float_range(0.0, 1.0),
            help="Normal VAF Cutoff in decimal format. Only sites BELOW this cutoff in normal will be considered.",
            default=0.02
        )
        self.parser.add_argument(
            '--tdna-vaf', type=float_range(0.0, 1.0),
            help="Tumor DNA VAF Cutoff in decimal format. Only sites above this cutoff will be considered.",
            default=0.25
        )
        self.parser.add_argument(
            '--trna-vaf', type=float_range(0.0, 1.0),
            help="Tumor RNA VAF Cutoff in decimal format. Only sites above this cutoff will be considered.",
            default=0.25
        )
        self.parser.add_argument(
            "--tumor-purity", type=float_range(0.0, 1.0),
            help="Value between 0 and 1 indicating the fraction of tumor cells in the tumor sample. Information is used during aggregate report creation for a simple estimation of whether variants are subclonal or clonal based on VAF. If not provided, purity is estimated directly from the VAFs.",
        )
        self.parser.add_argument(
            "--maximum-transcript-support-level", type=int,
            help="The threshold to use for filtering epitopes on the Ensembl transcript support level (TSL). "
                 + "Keep all epitopes with a transcript support level <= to this cutoff.",
            default=1,
            choices=[1, 2, 3, 4, 5]
        )
        self.parser.add_argument(
            "--biotypes", type=lambda s:[a for a in s.split(',')],
            help="A list of biotypes to use for pre-filtering transcripts for processing in the pipeline.",
            default=['protein_coding']
        )

    def genes_of_interest_args(self):
        self.parser.add_argument(
            '--genes-of-interest-file',
            help="A genes of interest file. Predictions resulting from variants on genes in this list will be marked in the result files. "
                 + "The file should be formatted to have each gene on a separate line without a header line. "
                 + "If no file is specified, the Cancer Gene Census list of high-confidence genes is used as the default."
        )

    def aggregated_report_args(self):
        self.parser.add_argument(
            '--aggregate-inclusion-binding-threshold', type=int,
            help="Threshold for including epitopes when creating the aggregate report",
            default=5000,
        )
        self.parser.add_argument(
            '--aggregate-inclusion-count-limit', type=int,
            help="Limit neoantigen candidates included in the aggregate report to only the best n candidates per variant.",
            default=15,
        )

    def pvacfuse(self):
        self.parser.add_argument(
            '--starfusion-file',
            help="Path to a star-fusion.fusion_predictions.tsv or star-fusion.fusion_predictions.abridged.tsv to extract read support and expression information from. Only used when running with AGFusion data."
        )
        self.parser.add_argument(
            '--read-support', type=int,
            help="Read Support Cutoff. Sites above this cutoff will be considered.",
            default=5
        )
        self.parser.add_argument(
            '--expn-val', type=float,
            help="Expression Cutoff. Expression is meassured as FFPM (fusion fragments per million total reads). Sites above this cutoff will be considered.",
            default=0.1
        )

    def pvacseq(self):
        self.parser.add_argument(
            "-p", "--phased-proximal-variants-vcf",
            help="A VCF with phased proximal variant information. Must be gzipped and tabix indexed."
        )
        self.parser.add_argument(
            "-c", "--minimum-fold-change", type=float,
            default=0.0,
            help="Minimum fold change between mutant (MT) binding score and wild-type (WT) score (fold change = WT/MT). "
                 + "The default is 0, which filters no results, but 1 is often a sensible choice "
                 + "(requiring that binding is better to the MT than WT peptide). "
                 + "This fold change is sometimes referred to as a differential agretopicity index.",
        )
        self.parser.add_argument(
            "--allele-specific-anchors",
            help="Use allele-specific anchor positions when tiering epitopes in the aggregate report. This option "
                 + "is available for 8, 9, 10, and 11mers and only for HLA-A, B, and C alleles. If this option is "
                 + "not enabled or as a fallback for unsupported lengths and alleles, the default positions of 1, "
                 + "2, epitope length - 1, and epitope length are used. Please see https://doi.org/10.1101/2020.12.08.416271 "
                 + "for more details.",
            default=False,
            action='store_true',
        )
        self.parser.add_argument(
            "--anchor-contribution-threshold", type=float_range(0.5,0.9),
            help="For determining allele-specific anchors, each position is assigned a score based on how binding is "
                 + "influenced by mutations. From these scores, the relative contribution of each position to the "
                 + "overall binding is calculated. Starting with the highest relative contribution, positions whose "
                 + "scores together account for the selected contribution threshold are assigned as anchor locations. "
                 + " As a result, a higher threshold leads to the inclusion of more positions to be considered anchors.",
            default=0.8
        )
        self.parser.add_argument(
            '--expn-val', type=float,
            default=1.0,
            help="Gene and Transcript Expression cutoff. Only sites above this cutoff will be considered.",
        )

    def pvacsplice(self):
        self.parser.add_argument(
            "annotated_vcf",
            help="A VEP-annotated single- or multi-sample VCF containing genotype and transcript information."
            + "The VCF may be gzipped (requires tabix index)."
        )
        self.parser.add_argument(
            "ref_fasta",
            help="A reference DNA FASTA file. Note: this input should be the same as the RegTools fasta input."
        )
        self.parser.add_argument(
            "gtf_file",
            help="A reference GTF file. Note: this input should be the same as the RegTools gtf input."
        )
        self.parser.add_argument(
            "-j", "--junction-score", type=int,
            help="Junction Coverage Cutoff. Only sites above this read depth cutoff will be considered.",
            default=10
        )
        self.parser.add_argument(
            "-v", "--variant-distance", type=int,
            help="Regulatory variants can lie inside or outside of splicing junction."
            + "Maximum distance window (upstream and downstream) for a variant outside the junction.",
            default=100
        )
        self.parser.add_argument(
            "-g", "--save-gtf",
            help="Save a tsv file from the uploaded filtered GTF data."
            + "Use this option to bypass GTF data upload time for multiple pVACsplice runs.",
            default=False,
            action='store_true'
        )
        self.parser.add_argument(
            "--anchor-types", type=pvacsplice_anchors(),
            help="The anchor types of junctions to use. Multiple anchors can be specified using a comma-separated list."
            + "Choices: A, D, NDA, DA, N",
            default=['A', 'D', 'NDA'],
        )
        # pvacsplice - filter on gene expression only (but keep txpn value in output)
        self.parser.add_argument(
            '--expn-val', type=float,
            default=1.0,
            help="Gene Expression cutoff. Only sites above this cutoff will be considered.",
        )

    def pvacvector(self):
        self.parser.add_argument(
            '-v', "--input-vcf",
            help="Path to original pVACseq input VCF file. Required if input file is a pVACseq TSV."
        )
        self.parser.add_argument(
            '-n', "--input-n-mer", default='25',
            help="Length of the peptide sequence to use when creating the FASTA from the pVACseq TSV.",
        )
        self.parser.add_argument(
            '--spacers', type=lambda s:[spacer for spacer in s.split(',')],
            help="Comma-separated list of spacers to use for testing junction epitopes. Include None to test junctions without spacers. Peptide combinations will be tested with each spacer in the order specified.",
            default="None,AAY,HHHH,GGS,GPGPG,HHAA,AAL,HH,HHC,HHH,HHHD,HHL,HHHC"
        )
        self.parser.add_argument(
            '--max-clip-length', type=int,
            help="Number of amino acids to permit clipping from the start and/or end of peptides in order to test novel junction epitopes when the first pass on the full peptide fails.",
            default=3,
        )
        self.parser.add_argument(
            '--allow-n-peptide-exclusion', type=int,
            help="If no solution is found after adding spacers and clipping peptides, attempt to find partial solutions with up to n peptides removed.",
            default=2,
        )

class PvacbindRunArgumentParser(RunArgumentParser):
    def __init__(self):
        tool_name = "pvacbind"
        input_file_help = "A FASTA file"
        RunArgumentParser.__init__(self, tool_name, input_file_help)
        self.epitope_args()
        self.binding_args(tool_name)
        self.prediction_args()
        self.aggregated_report_args()

class PvacfuseRunArgumentParser(RunArgumentParser):
    def __init__(self):
        tool_name = "pvacfuse"
        input_file_help="An AGFusion output directory or Arriba fusion.tsv output file."
        RunArgumentParser.__init__(self, tool_name, input_file_help)
        self.epitope_args()
        self.binding_args(tool_name)
        self.prediction_args()
        self.fasta_generation()
        self.genes_of_interest_args()
        self.aggregated_report_args()
        self.pvacfuse()

class PvacspliceRunArgumentParser(RunArgumentParser):
    def __init__(self):
        tool_name = "pvacsplice"
        input_file_help = "RegTools junctions output TSV file"
        RunArgumentParser.__init__(self, tool_name, input_file_help)
        self.epitope_args()
        self.binding_args(tool_name)
        self.pass_only_args()
        self.expression_coverage_args()
        self.prediction_args()
        self.genes_of_interest_args()
        self.aggregated_report_args()
        self.pvacsplice()


class PvacseqRunArgumentParser(RunArgumentParser):
    def __init__(self):
        tool_name = "pvacseq"
        input_file_help = (
            "A VEP-annotated single- or multi-sample VCF containing genotype, transcript, "
            "Wildtype protein sequence, and Frameshift protein sequence information."
            "The VCF may be gzipped (requires tabix index)."
        )
        RunArgumentParser.__init__(self, tool_name, input_file_help)
        self.epitope_args()
        self.binding_args(tool_name)
        self.pass_only_args()
        self.expression_coverage_args()
        self.prediction_args()
        self.fasta_generation()
        self.genes_of_interest_args()
        self.aggregated_report_args()
        self.pvacseq()

class PvacvectorRunArgumentParser(RunArgumentParser):
    def __init__(self):
        tool_name = 'pvacvector'
        input_file_help = "A .fa file with peptides or a pVACseq .tsv file with epitopes to use for vector design."
        RunArgumentParser.__init__(self, tool_name, input_file_help)
        self.parser.add_argument(
            "-e1", "--class-i-epitope-length", type=lambda s:[int(epl) for epl in s.split(',')],
            default=[8,9,10,11],
            help="Length of MHC Class I junctional epitopes to predict. "
                 + "Multiple epitope lengths can be specified using a comma-separated list. "
                 + "Typical epitope lengths vary between 8-15. "
                 + "Required for Class I prediction algorithms.",
        )
        self.parser.add_argument(
                "-e2", "--class-ii-epitope-length", type=lambda s:[int(epl) for epl in s.split(',')],
            default=[12,13,14,15,16,17,18],
            help="Length of MHC Class II junctional epitopes to predict. "
                 + "Multiple epitope lengths can be specified using a comma-separated list. "
                 + "Typical epitope lengths vary between 11-30. "
                 + "Required for Class II prediction algorithms.",
        )
        self.parser.add_argument(
            "-b","--binding-threshold", type=int,
            default=500,
            help="Fail junctions where any junctional epitope has ic50 binding scores below this value.",
        )
        self.parser.add_argument(
            '--percentile-threshold', type=float_range(0.0,100.0),
            help="Fail junctions where any junctional epitope "
                 +"has a percentile rank below this value."
        )
        self.parser.add_argument(
            '--percentile-threshold-strategy',
            choices=['conservative', 'exploratory'],
            help="Specify the how to evaluate junctional epitopes if a percentile threshold is set. "
                 + " The 'conservative' option fails a junction if a junctional epitope fails EITHER the binding threshold OR the percentile threshold (default)."
                 + " The 'exploratory' option fails a junction only if a junctional epitope fails BOTH the binding threshold AND the percentile threshold.",
            default="conservative",
        )
        self.parser.add_argument(
            '--allele-specific-binding-thresholds',
            help="Use allele-specific binding thresholds when evaluating junctional epitopes. To print the allele-specific binding thresholds run `pvacvector allele_specific_cutoffs`. "
                 + "If an allele does not have a special threshold value, the `--binding-threshold` value will be used.",
            default=False,
            action='store_true',
        )
        self.parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            default='median',
            help="The ic50 scoring metric to use when evaluating junctional epitopes by binding-threshold. "
                 + "lowest: Use the best MT Score (i.e. the lowest MT ic50 binding score of all chosen prediction methods). "
                 + "median: Use the median MT Score (i.e. the  median MT ic50 binding score of all chosen prediction methods)."
        )
        self.pvacvector()
