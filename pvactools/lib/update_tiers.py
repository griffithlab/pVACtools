from abc import ABCMeta, abstractmethod
import os
import csv
import ast
import pandas as pd
import tempfile
import shutil

from pvactools.lib.prediction_class import PredictionClass
from pvactools.lib.run_utils import get_anchor_positions, is_preferred_transcript

class UpdateTiers:
    def __init__(self):
        self.hla_types = pd.read_csv(self.input_file, delimiter="\t", usecols=["Allele"])['Allele'].unique()
        allele_specific_binding_thresholds = {}
        for hla_type in self.hla_types:
            threshold = PredictionClass.cutoff_for_allele(hla_type)
            if threshold is not None:
                allele_specific_binding_thresholds[hla_type] = float(threshold)
        self.allele_specific_binding_thresholds = allele_specific_binding_thresholds

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
            aggregate_inclusion_binding_threshold=5000,
            aggregate_inclusion_count_limit=15,
        ):
        self.input_file = input_file
        self.output_file = tempfile.NamedTemporaryFile()
        self.vaf_clonal = vaf_clonal
        self.binding_threshold = binding_threshold
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
        self.percentile_threshold = percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.aggregate_inclusion_binding_threshold = aggregate_inclusion_binding_threshold
        self.aggregate_inclusion_count_limit = aggregate_inclusion_count_limit
        self.allele_expr_threshold = trna_vaf * expn_val * 10
        self.trna_cov = trna_cov
        self.trna_vaf = trna_vaf
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level
        anchor_probabilities = {}
        for length in [8, 9, 10, 11]:
            base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
            file_name = os.path.join(base_dir, 'tools', 'pvacview', 'data', "Normalized_anchor_predictions_{}_mer.tsv".format(length))
            probs = {}
            with open(file_name, 'r') as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for line in reader:
                    hla = line.pop('HLA')
                    probs[hla] = line
            anchor_probabilities[length] = probs
        self.anchor_probabilities = anchor_probabilities

        mouse_anchor_positions = {}
        for length in [8, 9, 10, 11]:
            base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
            file_name = os.path.join(base_dir, 'tools', 'pvacview', 'data', "mouse_anchor_predictions_{}_mer.tsv".format(length))
            values = {}
            with open(file_name, 'r') as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for line in reader:
                    allele = line.pop('Allele')
                    values[allele] = {int(k): ast.literal_eval(v) for k, v in line.items()}
            mouse_anchor_positions[length] = values
        self.mouse_anchor_positions = mouse_anchor_positions

        self.allele_specific_anchors = allele_specific_anchors
        self.anchor_contribution_threshold = anchor_contribution_threshold
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

        anchor_residue_pass = self.is_anchor_residue_pass(mutation, binding_threshold)

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

    def is_anchor_residue_pass(self, mutation, binding_threshold):
        anchors = get_anchor_positions(mutation['Allele'], len(mutation['Best Peptide']), self.allele_specific_anchors, self.anchor_probabilities, self.anchor_contribution_threshold, self.mouse_anchor_positions)
        # parse out mutation positions from str
        position = mutation["Pos"]
        if position == 'NA':
            return True
        else:
            positions = position.split(", ")
            if len(positions) > 2:
                return True
            anchor_residue_pass = True
            if all(int(pos) in anchors for pos in positions):
                if mutation["IC50 WT"] == 'NA':
                    anchor_residue_pass = False
                elif float(mutation["IC50 WT"]) < binding_threshold:
                    anchor_residue_pass = False
            return anchor_residue_pass

    def sort_table(self, output_lines):
        #make sure the tiers sort in the expected order
        df = pd.DataFrame.from_dict(output_lines)
        tier_sorter = ["Pass", "PoorBinder", "RefMatch", "PoorTranscript", "LowExpr", "Anchor", "Subclonal", "ProbPos", "Poor", "NoExpr"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        df["rank_tier"] = df['Tier'].map(sorter_index)

        df["rank_ic50"] = pd.to_numeric(df["IC50 MT"]).rank(ascending=True, method='dense')
        df["rank_expr"] = pd.to_numeric(df["Allele Expr"], errors='coerce').rank(ascending=False, method='dense', na_option="bottom")
        df["rank"] = df["rank_ic50"] + df["rank_expr"]

        df.sort_values(by=["rank_tier", "rank", "Gene", "AA Change"], inplace=True, ascending=True)

        df.drop(labels='rank_tier', axis=1, inplace=True)
        df.drop(labels='rank_ic50', axis=1, inplace=True)
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
    ):
        self.input_file = input_file
        self.output_file = tempfile.NamedTemporaryFile()
        self.binding_threshold = binding_threshold
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
        self.percentile_threshold = percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.read_support = read_support
        self.expn_val = expn_val
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

        df['ic50_num'] = pd.to_numeric(df['IC50 MT'])
        df["rank_ic50"] = df["ic50_num"].rank(ascending=True, method='dense')
        df["rank_expr"] = pd.to_numeric(df["Expr"], errors='coerce').rank(ascending=False, method='dense', na_option="bottom")
        df["rank"] = df["rank_ic50"] + df["rank_expr"]

        df.sort_values(by=["rank_tier", "rank", "ic50_num", "ID"], inplace=True, ascending=True)

        df.drop(labels='rank_tier', axis=1, inplace=True)
        df.drop(labels='ic50_num', axis=1, inplace=True)
        df.drop(labels='rank_ic50', axis=1, inplace=True)
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

        df["rank_ic50"] = pd.to_numeric(df["IC50 MT"]).rank(ascending=True, method='dense')
        df["rank_expr"] = pd.to_numeric(df["Allele Expr"], errors='coerce').rank(ascending=False, method='dense', na_option="bottom")
        df["rank"] = df["rank_ic50"] + df["rank_expr"]

        df.sort_values(by=["rank_tier", "rank", "Gene", "Transcript", "AA Change"], inplace=True, ascending=True)

        df.drop(labels='rank_tier', axis=1, inplace=True)
        df.drop(labels='rank_ic50', axis=1, inplace=True)
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
            top_score_metric="median",
        ):
        self.input_file = input_file
        self.output_file = tempfile.NamedTemporaryFile()
        self.binding_threshold = binding_threshold
        self.percentile_threshold = percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
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

        df['ic50_num'] = pd.to_numeric(df['IC50 MT'])

        df.sort_values(by=["rank_tier", "ic50_num", "ID"], inplace=True, ascending=[True, True, True])

        df.drop(labels='rank_tier', axis=1, inplace=True)
        df.drop(labels='ic50_num', axis=1, inplace=True)
        return df

