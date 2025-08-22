from abc import ABCMeta, abstractmethod
import os
import csv
import ast
import pandas as pd
import tempfile
import shutil
import argparse

from pvactools.lib.prediction_class import PredictionClass
from pvactools.lib.run_utils import is_preferred_transcript, float_range, transcript_prioritization_strategy
from pvactools.lib.anchor_residue_pass import AnchorResiduePass

class UpdateTiers:
    def __init__(self):
        self.hla_types = pd.read_csv(self.input_file, delimiter="\t", usecols=["Allele"])['Allele'].unique()
        allele_specific_binding_thresholds = {}
        for hla_type in self.hla_types:
            threshold = PredictionClass.cutoff_for_allele(hla_type)
            if threshold is not None:
                allele_specific_binding_thresholds[hla_type] = float(threshold)
        self.allele_specific_binding_thresholds = allele_specific_binding_thresholds
        super().__init__()

    def execute(self):
        with open(self.input_file) as input_fh:
            reader = csv.DictReader(input_fh, delimiter="\t")
            output_lines = []
            for line in reader:
                line['Tier'] = self.get_tier(line)
                output_lines.append(line)
        output_df = self.sort_table(output_lines)
        output_df.to_csv(self.output_file, sep='\t', na_rep='NA', index=False, float_format='%.3f')
        shutil.copy(self.output_file.name, self.input_file)
        self.output_file.close()

    @abstractmethod
    def get_tier(self, mutation):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def sort_table(self, output_lines):
        raise Exception("Must implement method in child class")

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s update_tiers' % tool,
            description="Update tiers in an aggregated report in order to, for example, use different thresholds or account for problematic position or reference match information if run after initial pipeline run.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="Input aggregated file with tiers to update. This file will be overwritten with the output."
        )
        if tool in ['pvacseq', 'pvacsplice']:
            parser.add_argument(
                'vaf_clonal', type=float_range(0.0, 1.0),
                help="The RNA VAF threshold to determine whether a candidate is considered clonal. Any candidates with RNA VAF < vaf_clonal/2 will be considered subclonal."
            )
        parser.add_argument(
            "-b","--binding-threshold", type=int,
            default=500,
            help="IC50 Binding Threshold to consider when evaluting the binding criteria. Candidates  where the mutant allele has ic50 binding scores below this value will be considered good binders.",
        )
        parser.add_argument(
            '--allele-specific-binding-thresholds',
            help="Use allele-specific binding thresholds when evaluatng the binding criteria for tiering. To print the allele-specific binding thresholds run `%s allele_specific_cutoffs`. " % tool
                 + "If an allele does not have a special threshold value, the `--binding-threshold` value will be used.",
            default=False,
            action='store_true',
        )
        parser.add_argument(
            '--percentile-threshold', type=float_range(0.0,100.0),
            help="Account for the IC50 percentile rank when evaluating the binding criteria for tiering. A candidate's "
                 +"percentile rank must be below this value."
        )
        parser.add_argument(
            '--percentile-threshold-strategy',
            choices=['conservative', 'exploratory'],
            help="Specify the candidate inclusion strategy. The 'conservative' option requires a candidate to pass BOTH the binding threshold and percentile threshold (default) in order to pass the binding criteria."
                 + " The 'exploratory' option requires a candidate to pass EITHER the binding threshold or the percentile threshold.",
            default="conservative",
        )
        parser.add_argument(
            '-m2', '--top-score-metric2',
            choices=['ic50','percentile'],
            default='ic50',
            help='Whether to use IC50 MT or to use %%ile MT score column when sorting candidates within a tier.'
        )
        if tool in ['pvacseq', 'pvacsplice']:
            parser.add_argument(
                '--trna-vaf', type=float_range(0.0, 1.0),
                help="Tumor RNA VAF Cutoff in decimal format to consider when evaluating the expression criteria. Only sites above this cutoff will be considered.",
                default=0.25
            )
            parser.add_argument(
                '--trna-cov', type=int,
                help="Tumor RNA Coverage Cutoff to consider when evaluating the expression criteria. Only sites above this read depth cutoff will be considered.",
                default=10
            )
        if tool in ['pvacseq', 'pvacfuse', 'pvacsplice']:
            parser.add_argument(
                '--expn-val', type=float,
                help="Expression Cutoff to consider when evaluating the expression criteria. Expression is meassured as FFPM (fusion fragments per million total reads). Sites above this cutoff will be considered.",
                default=0.1
            )
        if tool == 'pvacfuse':
            parser.add_argument(
                '--read-support', type=int,
                help="Read Support Cutoff. Sites above this cutoff will be considered.",
                default=5
            )
        if tool in ['pvacseq', 'pvacsplice']:
            parser.add_argument(
                "--transcript-prioritization-strategy", type=transcript_prioritization_strategy(),
                help="Specify the criteria to consider when evaluating transcripts of the neoantigen candidates. "
                     + "'canonical' will consider a candidate to come from a good transcript if the transcript is a Ensembl canonical transcript. "
                     + "'mane_select' will consider a candidate to come from a good transcript if the transcript is a MANE select transcript. "
                     + "'tsl' will consider a candidate to come from a good transcript if the transcript's support level (TSL) passes the --maximum-transcript-support-level. "
                     + "When selecting more than one criteria, a transcript meeting EITHER of the selected criteria will be prioritized/selected.",
                default=['canonical', 'mane_select', 'tsl']
            )
            parser.add_argument(
                "--maximum-transcript-support-level", type=int,
                help="The threshold to use for filtering epitopes on the Ensembl transcript support level (TSL). "
                     + "Keep all epitopes with a transcript support level <= to this cutoff.",
                default=1,
                choices=[1, 2, 3, 4, 5]
            )
        if tool == 'pvacseq':
            parser.add_argument(
                "--allele-specific-anchors",
                help="Use allele-specific anchor positions when evaluating the anchor criteria for tiering epitopes in the aggregate report. This option "
                     + "is available for 8, 9, 10, and 11mers and only for HLA-A, B, and C alleles. If this option is "
                     + "not enabled or as a fallback for unsupported lengths and alleles, the default positions of 1, "
                     + "2, epitope length - 1, and epitope length are used. Please see https://doi.org/10.1101/2020.12.08.416271 "
                     + "for more details.",
                default=False,
                action='store_true',
            )
            parser.add_argument(
                "--anchor-contribution-threshold", type=float_range(0.5,0.9),
                help="For determining allele-specific anchors, each position is assigned a score based on how binding is "
                     + "influenced by mutations. From these scores, the relative contribution of each position to the "
                     + "overall binding is calculated. Starting with the highest relative contribution, positions whose "
                     + "scores together account for the selected contribution threshold are assigned as anchor locations. "
                     + " As a result, a higher threshold leads to the inclusion of more positions to be considered anchors.",
                default=0.8
            )
        return parser

class PvacseqUpdateTiers(UpdateTiers, metaclass=ABCMeta):
    def __init__(
            self,
            input_file,
            vaf_clonal,
            binding_threshold=500,
            trna_vaf=0.25,
            trna_cov=10,
            expn_val=1,
            transcript_prioritization_strategy=['mane_select', 'canonical', 'tsl'],
            maximum_transcript_support_level=1,
            percentile_threshold=None,
            percentile_threshold_strategy='conservative',
            allele_specific_binding_thresholds=False,
            allele_specific_anchors=False,
            anchor_contribution_threshold=0.8,
            top_score_metric2='ic50',
        ):
        self.input_file = input_file
        self.output_file = tempfile.NamedTemporaryFile()
        self.vaf_clonal = vaf_clonal
        self.binding_threshold = binding_threshold
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
        self.percentile_threshold = percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.allele_expr_threshold = trna_vaf * expn_val * 10
        self.trna_cov = trna_cov
        self.trna_vaf = trna_vaf
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level
        if top_score_metric2 == "percentile":
            self.top_score_mode = "%ile MT"
        else:
            self.top_score_mode = "IC50 MT"
        super().__init__()
        self.anchor_calculator = AnchorResiduePass(binding_threshold, self.use_allele_specific_binding_thresholds, self.allele_specific_binding_thresholds, allele_specific_anchors, anchor_contribution_threshold)

    def get_tier(self, mutation):
        if self.use_allele_specific_binding_thresholds and mutation['Allele'] in self.allele_specific_binding_thresholds:
            binding_threshold = self.allele_specific_binding_thresholds[mutation['Allele']]
        else:
            binding_threshold = self.binding_threshold

        ic50_pass = float(mutation["IC50 MT"]) < binding_threshold
        percentile_pass = (
            self.percentile_threshold is None or
            float(mutation["%ile MT"]) < self.percentile_threshold
        )
        binding_pass = (
            (ic50_pass and percentile_pass)
            if self.percentile_threshold_strategy == 'conservative'
            else (ic50_pass or percentile_pass)
        )

        anchor_residue_pass = self.anchor_calculator.is_anchor_residue_pass(mutation)

        transcript_pass = is_preferred_transcript(mutation, self.transcript_prioritization_strategy, self.maximum_transcript_support_level)

        refmatch_pass = True
        if 'Ref Match' in mutation:
            refmatch_pass = mutation['Ref Match'] != "True"

        probaa_pass = True
        if 'Prob Pos' in mutation:
            probaa_pass = mutation['Prob Pos'] == 'None'

        allele_expr_pass = True
        if (mutation['Allele Expr'] != 'NA' and float(mutation['Allele Expr']) <= self.allele_expr_threshold):
            allele_expr_pass = False

        vaf_clonal_pass = True
        if (mutation['DNA VAF'] != 'NA' and float(mutation['DNA VAF']) < (self.vaf_clonal/2)):
            vaf_clonal_pass = False

        #writing these out as explicitly as possible for ease of understanding
        if (binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           transcript_pass and
           anchor_residue_pass and
           refmatch_pass and
           probaa_pass):
            return "Pass"

        #poor binder
        if (not binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           transcript_pass and
           anchor_residue_pass and
           refmatch_pass and
           probaa_pass):
            return "PoorBinder"

        #has reference match
        if (binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           transcript_pass and
           anchor_residue_pass and
           not refmatch_pass and
           probaa_pass):
            return "RefMatch"

        #has problematic positions
        if (binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           transcript_pass and
           anchor_residue_pass and
           refmatch_pass and
           not probaa_pass):
            return "ProbPos"

        #transcript doesn't match the prioritization criteria
        if (binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           not transcript_pass and
           refmatch_pass and
           probaa_pass):
            return "PoorTranscript"

        #anchor residues
        if (binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           transcript_pass and
           not anchor_residue_pass and
           refmatch_pass and
           probaa_pass):
            return "Anchor"

        #not in founding clone
        if (binding_pass and
           allele_expr_pass and
           not vaf_clonal_pass and
           transcript_pass and
           anchor_residue_pass and
           refmatch_pass and
           probaa_pass):
            return "Subclonal"

        #relax expression.  Include sites that have reasonable vaf but zero overall gene expression
        lowexpr=False
        if mutation['RNA VAF'] != 'NA' and mutation['RNA Expr'] != 'NA' and ['RNA Depth'] != 'NA' and mutation['Allele Expr'] != 'NA':
            if ((float(mutation["Allele Expr"]) > 0) or
               (float(mutation["RNA Expr"]) == 0 and
               float(mutation["RNA Depth"]) > self.trna_cov and
               float(mutation["RNA VAF"]) > self.trna_vaf)):
                lowexpr=True

        #if low expression is the only strike against it, it gets lowexpr label (multiple strikes will pass through to poor)
        if (binding_pass and
           lowexpr and
           vaf_clonal_pass and
           transcript_pass and
           anchor_residue_pass and
           refmatch_pass and
           probaa_pass):
            return "LowExpr"

        #zero expression
        if (((mutation['RNA Expr'] != 'NA' and float(mutation["RNA Expr"]) == 0) or
            (mutation['RNA VAF'] != 'NA' and float(mutation["RNA VAF"]) == 0)) and
            not lowexpr):
            return "NoExpr"

        #everything else
        return "Poor"

    def sort_table(self, output_lines):
        #make sure the tiers sort in the expected order
        df = pd.DataFrame.from_dict(output_lines)
        tier_sorter = ["Pass", "PoorBinder", "RefMatch", "PoorTranscript", "LowExpr", "Anchor", "Subclonal", "ProbPos", "Poor", "NoExpr"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        df["rank_tier"] = df['Tier'].map(sorter_index)

        df["rank_binding"] = pd.to_numeric(df[self.top_score_mode]).rank(ascending=True, method='dense')
        df["rank_expr"] = pd.to_numeric(df["Allele Expr"], errors='coerce').rank(ascending=False, method='dense', na_option="bottom")
        df["rank"] = df["rank_binding"] + df["rank_expr"]

        df.sort_values(by=["rank_tier", "rank", "Gene", "AA Change"], inplace=True, ascending=True)

        df.drop(labels='rank_tier', axis=1, inplace=True)
        df.drop(labels='rank_binding', axis=1, inplace=True)
        df.drop(labels='rank_expr', axis=1, inplace=True)
        df.drop(labels='rank', axis=1, inplace=True)

        return df

class PvacfuseUpdateTiers(UpdateTiers, metaclass=ABCMeta):
    def __init__(
        self,
        input_file,
        binding_threshold=500,
        percentile_threshold=None,
        percentile_threshold_strategy='conservative',
        allele_specific_binding_thresholds=False,
        read_support=5,
        expn_val=0.1,
        top_score_metric2="ic50"
    ):
        self.input_file = input_file
        self.output_file = tempfile.NamedTemporaryFile()
        self.binding_threshold = binding_threshold
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
        self.percentile_threshold = percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.read_support = read_support
        self.expn_val = expn_val
        if top_score_metric2 == "percentile":
            self.top_score_mode = "%ile MT"
        else:
            self.top_score_mode = "IC50 MT"
        super().__init__()

    def get_tier(self, mutation):
        if self.use_allele_specific_binding_thresholds and mutation['Allele'] in self.allele_specific_binding_thresholds:
            binding_threshold = self.allele_specific_binding_thresholds[mutation['Allele']]
        else:
            binding_threshold = self.binding_threshold

        ic50_pass = float(mutation["IC50 MT"]) < binding_threshold
        percentile_pass = (
            self.percentile_threshold is None or
            float(mutation["%ile MT"]) < self.percentile_threshold
        )
        binding_pass = (
            (ic50_pass and percentile_pass)
            if self.percentile_threshold_strategy == 'conservative'
            else (ic50_pass or percentile_pass)
        )

        low_read_support = False
        if mutation['Read Support'] != 'NA' and float(mutation['Read Support']) < self.read_support:
            low_read_support = True

        low_expr = False
        if mutation['Expr'] != 'NA' and float(mutation['Expr']) < self.expn_val:
            low_expr = True

        refmatch_pass = True
        if 'Ref Match' in mutation:
            refmatch_pass = mutation['Ref Match'] != "True"

        probaa_pass = True
        if 'Prob Pos' in mutation:
            probaa_pass = mutation['Prob Pos'] == 'None'

        if (binding_pass and
          not low_read_support and
          not low_expr and
          refmatch_pass and
          probaa_pass):
            return "Pass"

        #poor binder
        if (not binding_pass and
          not low_read_support and
          not low_expr and
          refmatch_pass and
          probaa_pass):
            return "PoorBinder"

        #has reference match
        if (binding_pass and
          not low_read_support and
          not low_expr and
          not refmatch_pass and
          probaa_pass):
            return "RefMatch"

        #has problematic positions
        if (binding_pass and
          not low_read_support and
          not low_expr and
          refmatch_pass and
          not probaa_pass):
            return "ProbPos"

        #low read support
        if (binding_pass and
          low_read_support and
          not low_expr and
          refmatch_pass and
          probaa_pass):
            return "LowReadSupport"

        #low expression
        if (binding_pass and
          not low_read_support and
          low_expr and
          refmatch_pass and
          probaa_pass):
            return "LowExpr"

        return "Poor"

    def sort_table(self, output_lines):
        df = pd.DataFrame.from_dict(output_lines)
        tier_sorter = ["Pass", "PoorBinder", "RefMatch", "LowReadSupport", "LowExpr", "ProbPos", "Poor"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        df["rank_tier"] = df['Tier'].map(sorter_index)

        df['binding_num'] = pd.to_numeric(df[self.top_score_mode])
        df["rank_binding"] = df["binding_num"].rank(ascending=True, method='dense')
        df["rank_expr"] = pd.to_numeric(df["Expr"], errors='coerce').rank(ascending=False, method='dense', na_option="bottom")
        df["rank"] = df["rank_binding"] + df["rank_expr"]

        df.sort_values(by=["rank_tier", "rank", "binding_num", "ID"], inplace=True, ascending=True)

        df.drop(labels='rank_tier', axis=1, inplace=True)
        df.drop(labels='binding_num', axis=1, inplace=True)
        df.drop(labels='rank_binding', axis=1, inplace=True)
        df.drop(labels='rank_expr', axis=1, inplace=True)
        df.drop(labels='rank', axis=1, inplace=True)
        return df

class PvacspliceUpdateTiers(UpdateTiers, metaclass=ABCMeta):
    def __init__(
        self,
        input_file,
        vaf_clonal,
        binding_threshold=500,
        allele_specific_binding_thresholds=False,
        percentile_threshold=None,
        percentile_threshold_strategy='conservative',
        trna_vaf=0.25,
        trna_cov=10,
        expn_val=1,
        transcript_prioritization_strategy=['mane_select', 'canonical', 'tsl'],
        maximum_transcript_support_level=1,
        top_score_metric2="ic50"
    ):
        self.input_file = input_file
        self.output_file = tempfile.NamedTemporaryFile()
        self.binding_threshold = binding_threshold
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
        self.percentile_threshold = percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.vaf_clonal = vaf_clonal
        self.allele_expr_threshold = trna_vaf * expn_val * 10
        self.trna_vaf = trna_vaf
        self.trna_cov = trna_cov
        self.expn_val = expn_val
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level
        if top_score_metric2 == "percentile":
            self.top_score_mode = "%ile MT"
        else:
            self.top_score_mode = "IC50 MT"
        super().__init__()

    def get_tier(self, mutation):
        if self.use_allele_specific_binding_thresholds and mutation['Allele'] in self.allele_specific_binding_thresholds:
            binding_threshold = self.allele_specific_binding_thresholds[mutation['Allele']]
        else:
            binding_threshold = self.binding_threshold

        ic50_pass = float(mutation["IC50 MT"]) < binding_threshold
        percentile_pass = (
            self.percentile_threshold is None or
            float(mutation["%ile MT"]) < self.percentile_threshold
        )
        binding_pass = (
            (ic50_pass and percentile_pass)
            if self.percentile_threshold_strategy == 'conservative'
            else (ic50_pass or percentile_pass)
        )

        transcript_pass = is_preferred_transcript(mutation, self.transcript_prioritization_strategy, self.maximum_transcript_support_level)

        refmatch_pass = True
        if 'Ref Match' in mutation:
            refmatch_pass = mutation['Ref Match'] != "True"

        probaa_pass = True
        if 'Prob Pos' in mutation:
            probaa_pass = mutation['Prob Pos'] == 'None'

        allele_expr_pass = True
        if (mutation['Allele Expr'] != 'NA' and float(mutation['Allele Expr']) <= self.allele_expr_threshold):
            allele_expr_pass = False

        vaf_clonal_pass = True
        if (mutation['DNA VAF'] != 'NA' and float(mutation['DNA VAF']) < (self.vaf_clonal/2)):
            vaf_clonal_pass = False

        #writing these out as explicitly as possible for ease of understanding
        if (binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           transcript_pass and
           refmatch_pass and
           probaa_pass):
            return "Pass"

        #poor binder
        if (not binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           transcript_pass and
           refmatch_pass and
           probaa_pass):
            return "PoorBinder"

        #has reference match
        if (binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           transcript_pass and
           not refmatch_pass and
           probaa_pass):
            return "RefMatch"

        #has problematic positions
        if (binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           transcript_pass and
           refmatch_pass and
           not probaa_pass):
            return "ProbPos"

        #transcript doesn't match the prioritization criteria
        if (binding_pass and
           allele_expr_pass and
           vaf_clonal_pass and
           not transcript_pass and
           refmatch_pass and
           probaa_pass):
            return "PoorTranscript"

        #not in founding clone
        if (binding_pass and
           allele_expr_pass and
           not vaf_clonal_pass and
           transcript_pass and
           refmatch_pass and
           probaa_pass):
            return "Subclonal"

        #relax expression.  Include sites that have reasonable vaf but zero overall gene expression
        lowexpr=False
        if mutation['RNA VAF'] != 'NA' and mutation['RNA Expr'] != 'NA' and ['RNA Depth'] != 'NA' and mutation['Allele Expr'] != 'NA':
            if ((float(mutation["Allele Expr"]) > 0) or
               (float(mutation["RNA Expr"]) == 0 and
               float(mutation["RNA Depth"]) > self.trna_cov and
               float(mutation["RNA VAF"]) > self.trna_vaf)):
                lowexpr=True

        #if low expression is the only strike against it, it gets lowexpr label (multiple strikes will pass through to poor)
        if (binding_pass and
           lowexpr and
           vaf_clonal_pass and
           transcript_pass and
           refmatch_pass and
           probaa_pass):
            return "LowExpr"

        #zero expression
        if (((mutation['RNA Expr'] != 'NA' and float(mutation["RNA Expr"]) == 0) or
            (mutation['RNA VAF'] != 'NA' and float(mutation["RNA VAF"]) == 0)) and
            not lowexpr):
            return "NoExpr"

        #everything else
        return "Poor"

    def sort_table(self, output_lines):
        df = pd.DataFrame.from_dict(output_lines)

        #make sure the tiers sort in the expected order
        tier_sorter = ["Pass", "PoorBinder", "RefMatch", "PoorTranscript", "LowExpr", "Subclonal", "ProbPos", "Poor", "NoExpr"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        df["rank_tier"] = df['Tier'].map(sorter_index)

        df["rank_binding"] = pd.to_numeric(df[self.top_score_mode]).rank(ascending=True, method='dense')
        df["rank_expr"] = pd.to_numeric(df["Allele Expr"], errors='coerce').rank(ascending=False, method='dense', na_option="bottom")
        df["rank"] = df["rank_binding"] + df["rank_expr"]

        df.sort_values(by=["rank_tier", "rank", "Gene", "Transcript", "AA Change"], inplace=True, ascending=True)

        df.drop(labels='rank_tier', axis=1, inplace=True)
        df.drop(labels='rank_binding', axis=1, inplace=True)
        df.drop(labels='rank_expr', axis=1, inplace=True)
        df.drop(labels='rank', axis=1, inplace=True)

        return df

class PvacbindUpdateTiers(UpdateTiers, metaclass=ABCMeta):
    def __init__(
            self,
            input_file,
            binding_threshold=500,
            percentile_threshold=None,
            percentile_threshold_strategy='conservative',
            allele_specific_binding_thresholds=False,
            top_score_metric2="ic50"
        ):
        self.input_file = input_file
        self.output_file = tempfile.NamedTemporaryFile()
        self.binding_threshold = binding_threshold
        self.percentile_threshold = percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
        if top_score_metric2 == "percentile":
            self.top_score_mode = "%ile MT"
        else:
            self.top_score_mode = "IC50 MT"
        super().__init__()

    def get_tier(self, mutation):
        if self.use_allele_specific_binding_thresholds and mutation['Allele'] in self.allele_specific_binding_thresholds:
            binding_threshold = self.allele_specific_binding_thresholds[mutation['Allele']]
        else:
            binding_threshold = self.binding_threshold

        ic50_pass = float(mutation["IC50 MT"]) < binding_threshold
        percentile_pass = (
            self.percentile_threshold is None or
            float(mutation["%ile MT"]) < self.percentile_threshold
        )
        binding_pass = (
            (ic50_pass and percentile_pass)
            if self.percentile_threshold_strategy == 'conservative'
            else (ic50_pass or percentile_pass)
        )

        refmatch_pass = True
        if 'Ref Match' in mutation:
            refmatch_pass = mutation['Ref Match'] != "True"

        probaa_pass = True
        if 'Prob Pos' in mutation:
            probaa_pass = mutation['Prob Pos'] == 'None'

        if (binding_pass and
            refmatch_pass and
            probaa_pass):
            return "Pass"

        #poor binder
        if (not binding_pass and
            refmatch_pass and
            probaa_pass):
            return "PoorBinder"

        #has reference match
        if (binding_pass and
           not refmatch_pass and
           probaa_pass):
            return "RefMatch"

        #has problematic positions
        if (binding_pass and
           refmatch_pass and
           not probaa_pass):
            return "ProbPos"

        return "Poor"

    def sort_table(self, output_lines):
        df = pd.DataFrame.from_dict(output_lines)

        tier_sorter = ["Pass", "PoorBinder", "RefMatch", "ProbPos", "Poor"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        df["rank_tier"] = df['Tier'].map(sorter_index)

        df['binding_num'] = pd.to_numeric(df[self.top_score_mode])

        df.sort_values(by=["rank_tier", "binding_num", "ID"], inplace=True, ascending=[True, True, True])

        df.drop(labels='rank_tier', axis=1, inplace=True)
        df.drop(labels='binding_num', axis=1, inplace=True)
        return df
