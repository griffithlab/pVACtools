from abc import ABCMeta
import argparse
from .prediction_class import *
import lib

class RunArgumentParser(metaclass=ABCMeta):
    def __init__(self, run_name, input_file_help):
        parser = argparse.ArgumentParser(run_name)

        parser.add_argument(
            "input_file",
            help=input_file_help
        )
        parser.add_argument(
            "sample_name",
            help="The name of the sample being processed. This will be used as a prefix for output files"
        )
        parser.add_argument(
            "allele", type=lambda s:[a for a in s.split(',')],
            help="Name of the allele to use for epitope prediction. "
                 + "Multiple alleles can be specified using a comma-separated list. "
                 + "For a list of available alleles, use: `pvacseq valid_alleles`",
        )
        parser.add_argument(
            "prediction_algorithms",
            choices=PredictionClass.prediction_methods(),
            nargs="+",
            help="The epitope prediction algorithms to use. Multiple prediction algorithms can be specified, separated by spaces",
        )
        parser.add_argument(
            "output_dir",
            help="The directory for writing all result files"
        )
        parser.add_argument(
            "-e", "--epitope-length", type=lambda s:[int(epl) for epl in s.split(',')],
            help="Length of subpeptides (neoepitopes) to predict. "
                 + "Multiple epitope lengths can be specified using a comma-separated list. "
                 + "Typical epitope lengths vary between 8-11. "
                 + "Required for Class I prediction algorithms",
        )
        parser.add_argument(
            "-l", "--peptide-sequence-length", type=int,
            default=21,
            help="Length of the peptide sequence to use when creating the FASTA. Default: 21",
        )
        parser.add_argument(
            "--iedb-install-directory",
            help="Directory that contains the local installation of IEDB MHC I and/or MHC II"
        )
        parser.add_argument(
            '--net-chop-method',
            choices=lib.net_chop.methods,
            default=None,
            help="NetChop prediction method to use (\"cterm\" for C term 3.0, \"20s\" for 20S 3.0). ",
        )
        parser.add_argument(
            '--netmhc-stab',
            action='store_true',
            help="Run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes"
        )
        parser.add_argument(
            '-t', '--top-result-per-mutation',
            action='store_true',
            help='Output only the top scoring result for each allele-peptide length combination for each variant. Default: False'
        )
        parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            default='median',
            help="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. "
                 + "lowest: Best MT Score/Corresponding Fold Change - lowest MT ic50 binding score/corresponding fold change of all chosen prediction methods. "
                 + "median: Median MT Score/Median Fold Change - median MT ic50 binding score/fold change of all chosen prediction methods. "
                 + "Default: median"
        )
        parser.add_argument(
            "-b","--binding-threshold", type=int,
            default=500,
            help="Report only epitopes where the mutant allele has ic50 binding scores below this value. Default: 500",
        )
        parser.add_argument(
            "-c", "--minimum-fold-change", type=int,
            default=0,
            help="Minimum fold change between mutant binding score and wild-type score. "
                 + "The default is 0, which filters no results, but 1 is often a sensible choice "
                 + "(requiring that binding is better to the MT than WT). Default: 0",
        )
        parser.add_argument(
            '--net-chop-threshold', type=float,
            default=0.5,
            help="NetChop prediction threshold. Default: 0.5",
        )
        parser.add_argument(
            '-a', '--additional-report-columns',
            choices=['sample_name'],
            help="Additional columns to output in the final report"
        )
        parser.add_argument(
            "-s", "--fasta-size",type=int,
            default=200,
            help="Number of fasta entries per IEDB request. "
                 + "For some resource-intensive prediction algorithms like Pickpocket and NetMHCpan it might be helpful to reduce this number. "
                 + "Needs to be an even number.",
        )
        parser.add_argument(
            "-r", "--iedb-retries",type=int,
            default=5,
            help="Number of retries when making requests to the IEDB RESTful web interface. Must be less than or equal to 100."
                 + "Default: 5"
        )
        parser.add_argument(
            "-d", "--downstream-sequence-length",
            default='1000',
            help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
                + "Use 'full' to include the full downstream sequence. Default: 1000"
        )
        parser.add_argument(
            "-k", "--keep-tmp-files",
            action='store_true',
            help="Keep intermediate output files. This migt be useful for debugging purposes.",
        )
        self.parser = parser

class PvacseqRunArgumentParser(RunArgumentParser):
    def __init__(self):
        run_name = "pvacseq run"
        input_file_help = (
            "A VEP-annotated single-sample VCF containing transcript, "
            "Wildtype protein sequence, and Downstream protein sequence information."
        )
        RunArgumentParser.__init__(self, run_name, input_file_help)

        self.parser.add_argument(
            "-i", "--additional-input-file-list",
            help="yaml file of additional files to be used as inputs, e.g. cufflinks output files. "
                 + "For an example of this yaml file run `pvacseq config_files additional_input_file_list`."
        )
        self.parser.add_argument(
            '--normal-cov', type=int,
            help="Normal Coverage Cutoff. Sites above this cutoff will be considered. " +
            "Default: 5",
            default=5
        )
        self.parser.add_argument(
            '--tdna-cov', type=int,
            help="Tumor DNA Coverage Cutoff. Sites above this cutoff will be considered. " +
            "Default: 10",
            default=10
        )
        self.parser.add_argument(
            '--trna-cov', type=int,
            help="Tumor RNA Coverage Cutoff. Sites above this cutoff will be considered. " +
            "Default: 10",
            default=10
        )
        self.parser.add_argument(
            '--normal-vaf', type=int,
            help="Normal VAF Cutoff. Sites BELOW this cutoff in normal will be considered. " +
            "Default: 2",
            default=2
        )
        self.parser.add_argument(
            '--tdna-vaf', type=int,
            help="Tumor DNA VAF Cutoff. Sites above this cutoff will be considered. " +
            "Default: 40",
            default=40
        )
        self.parser.add_argument(
            '--trna-vaf', type=int,
            help="Tumor RNA VAF Cutoff. Sites above this cutoff will be considered. " +
            "Default: 40",
            default=40
        )
        self.parser.add_argument(
            '--expn-val', type=int,
            default=1,
            help="Gene and Transcript Expression cutoff. Sites above this cutoff will be considered. Default: 1",
        )

class PvacfuseRunArgumentParser(RunArgumentParser):
    def __init__(self):
        run_name = "pvacfuse run"
        input_file_help = "A INTEGRATE-Neo bedpe file with fusions."
        RunArgumentParser.__init__(self, run_name, input_file_help)
