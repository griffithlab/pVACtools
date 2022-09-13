from abc import ABCMeta
import argparse

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
            "-e1", "--class-i-epitope-length", type=lambda s:[int(epl) for epl in s.split(',')],
            default=[8,9,10,11],
            help="Length of MHC Class I subpeptides (neoepitopes) to predict. "
                 + "Multiple epitope lengths can be specified using a comma-separated list. "
                 + "Typical epitope lengths vary between 8-15. "
                 + "Required for Class I prediction algorithms.",
        )
        parser.add_argument(
                "-e2", "--class-ii-epitope-length", type=lambda s:[int(epl) for epl in s.split(',')],
            default=[12,13,14,15,16,17,18],
            help="Length of MHC Class II subpeptides (neoepitopes) to predict. "
                 + "Multiple epitope lengths can be specified using a comma-separated list. "
                 + "Typical epitope lengths vary between 11-30. "
                 + "Required for Class II prediction algorithms.",
        )
        parser.add_argument(
            "--iedb-install-directory",
            help="Directory that contains the local installation of IEDB MHC I and/or MHC II."
        )
        parser.add_argument(
            "-b","--binding-threshold", type=int,
            default=500,
            help="Report only epitopes where the mutant allele has ic50 binding scores below this value.",
        )
        parser.add_argument(
            '--percentile-threshold', type=float,
            help="Report only epitopes where the mutant allele "
                 +"has a percentile rank below this value."
        )
        parser.add_argument(
            '--allele-specific-binding-thresholds',
            help="Use allele-specific binding thresholds. To print the allele-specific binding thresholds run `%s allele_specific_cutoffs`. " % tool_name
                 + "If an allele does not have a special threshold value, the `--binding-threshold` value will be used.",
            default=False,
            action='store_true',
        )
        parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            default='median',
            help="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. "
                 + "lowest: Use the best MT Score and Corresponding Fold Change (i.e. the lowest MT ic50 binding score and corresponding fold change of all chosen prediction methods). "
                 + "median: Use the median MT Score and Median Fold Change (i.e. the  median MT ic50 binding score and fold change of all chosen prediction methods)."
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
        self.parser = parser

class PredictionRunArgumentParser(RunArgumentParser):
    def __init__(self, tool_name, input_file_help):
        RunArgumentParser.__init__(self, tool_name, input_file_help)
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

class PvacbindRunArgumentParser(PredictionRunArgumentParser):
    def __init__(self):
        tool_name = "pvacbind"
        input_file_help = "A FASTA file"
        PredictionRunArgumentParser.__init__(self, tool_name, input_file_help)

class PredictionRunWithFastaGenerationArgumentParser(PredictionRunArgumentParser):
    def __init__(self, tool_name, input_file_help):
        PredictionRunArgumentParser.__init__(self, tool_name, input_file_help)
        self.parser.add_argument(
            "-d", "--downstream-sequence-length",
            default='1000',
            help="Cap to limit the downstream sequence length for frameshifts when creating the FASTA file. "
                + "Use 'full' to include the full downstream sequence."
        )

class PvacseqRunArgumentParser(PredictionRunWithFastaGenerationArgumentParser):
    def __init__(self):
        tool_name = "pvacseq"
        input_file_help = (
            "A VEP-annotated single- or multi-sample VCF containing genotype, transcript, "
            "Wildtype protein sequence, and Downstream protein sequence information."
            "The VCF may be gzipped (requires tabix index)."
        )
        PredictionRunWithFastaGenerationArgumentParser.__init__(self, tool_name, input_file_help)

        self.parser.add_argument(
            '--normal-sample-name',
            help="In a multi-sample VCF, the name of the matched normal sample."
        )
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
            '--normal-vaf', type=float_range(0.0,1.0),
            help="Normal VAF Cutoff in decimal format. Only sites BELOW this cutoff in normal will be considered.",
            default=0.02
        )
        self.parser.add_argument(
            '--tdna-vaf', type=float_range(0.0,1.0),
            help="Tumor DNA VAF Cutoff in decimal format. Only sites above this cutoff will be considered.",
            default=0.25
        )
        self.parser.add_argument(
            '--trna-vaf', type=float_range(0.0,1.0),
            help="Tumor RNA VAF Cutoff in decimal format. Only sites above this cutoff will be considered.",
            default=0.25
        )
        self.parser.add_argument(
            '--expn-val', type=float,
            default=1.0,
            help="Gene and Transcript Expression cutoff. Only sites above this cutoff will be considered.",
        )
        self.parser.add_argument(
            "--maximum-transcript-support-level", type=int,
            help="The threshold to use for filtering epitopes on the Ensembl transcript support level (TSL). "
            +"Keep all epitopes with a transcript support level <= to this cutoff.",
            default=1,
            choices=[1,2,3,4,5]
        )
        self.parser.add_argument(
            '--pass-only',
            help="Only process VCF entries with a PASS status.",
            default=False,
            action='store_true'
        )
        self.parser.add_argument(
            "--tumor-purity",
            help="Value between 0 and 1 indicating the fraction of tumor cells in the tumor sample. Information is used during aggregate report creation for a simple estimation of whether variants are subclonal or clonal based on VAF. If not provided, purity is estimated directly from the VAFs.",
            type=float,
        )


class PvacfuseRunArgumentParser(PredictionRunWithFastaGenerationArgumentParser):
    def __init__(self):
        tool_name = "pvacfuse"
        input_file_help = "An AGfusion output directory."
        PredictionRunWithFastaGenerationArgumentParser.__init__(self, tool_name, input_file_help)

class PvacvectorRunArgumentParser(RunArgumentParser):
    def __init__(self):
        tool_name = 'pvacvector'
        input_file_help = "A .fa file with peptides or a pVACseq .tsv file with epitopes to use for vector design."
        RunArgumentParser.__init__(self, tool_name, input_file_help)
        self.parser.add_argument(
            '-v', "--input_vcf",
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
