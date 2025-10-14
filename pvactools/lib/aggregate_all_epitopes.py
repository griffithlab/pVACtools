import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import json
import os
import shutil
from abc import ABCMeta, abstractmethod
import itertools
import csv
import glob
import ast
import math

from pvactools.lib.run_utils import is_preferred_transcript
from pvactools.lib.prediction_class import PredictionClass
from pvactools.lib.update_tiers import PvacseqUpdateTiers, PvacfuseUpdateTiers, PvacspliceUpdateTiers, PvacbindUpdateTiers
from pvactools.lib.anchor_residue_pass import AnchorResiduePass
from pvactools.lib.get_best_candidate import PvacseqBestCandidate, PvacfuseBestCandidate, PvacbindBestCandidate, PvacspliceBestCandidate

class AggregateAllEpitopes:
    def __init__(self):
        self.hla_types = pd.read_csv(self.input_file, delimiter="\t", usecols=["HLA Allele"])['HLA Allele'].unique()
        thresholds = {}
        for hla_type in self.hla_types:
            threshold = PredictionClass.cutoff_for_allele(hla_type)
            if threshold is not None:
                thresholds[hla_type] = float(threshold)
        self.allele_specific_binding_thresholds = thresholds

    @abstractmethod
    def get_list_unique_mutation_keys(self):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def calculate_clonal_vaf(self):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def read_input_file(self, used_columns, dtypes):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_sub_df(self, df, key):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_best_binder(self, df):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_included_df(self, df):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_good_binders_df(self, included_df):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_unique_peptide_hla_counts(self, included_df):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_included_df_metrics(self, included_df, prediction_algorithms, el_algorithms, percentile_algorithms):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def calculate_unique_peptide_count(self, included_df):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def calculate_good_binder_count(self, good_binders_df):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_default_annotation_count(self):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def assemble_result_line(self, best, key, hla, anno_count, included_peptide_count, good_binder_count):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_metrics(self, peptides, best):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def write_metrics_file(self, metrics):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def copy_pvacview_r_files(self):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def tier_aggregated_report(self):
        raise Exception("Must implement method in child class")

    def get_best_mut_line(self, df, key, prediction_algorithms, el_algorithms, percentile_algorithms):
        #order by best median score and get best ic50 peptide
        best = self.get_best_binder(df)

        #these are all lines meeting the aggregate inclusion binding threshold
        included_df = self.get_included_df(df)
        best_df = pd.DataFrame.from_dict([best])
        if not best_df.index.isin(included_df.index).all():
            included_df = pd.concat([included_df, best_df])
        best_df = best_df.to_dict()
        good_binders_df = self.get_good_binders_df(included_df)
        peptide_hla_counts = self.get_unique_peptide_hla_counts(good_binders_df)
        hla_counts = Counter(peptide_hla_counts["HLA Allele"])
        hla = dict(map(lambda x : (x, hla_counts[x]) if x in hla_counts else (x, ""), self.hla_types))

        #get a list of all unique gene/transcript/aa_change combinations
        #store a count of all unique peptides that passed
        (peptides, anno_count) = self.get_included_df_metrics(included_df, prediction_algorithms, el_algorithms, percentile_algorithms)
        included_peptide_count = self.calculate_unique_peptide_count(included_df)
        good_binder_count = self.calculate_good_binder_count(good_binders_df)

        #assemble the line
        out_dict = self.assemble_result_line(best, key, hla, anno_count, included_peptide_count, good_binder_count)

        metric = self.get_metrics(peptides, best)
        return (out_dict, metric)

    def determine_used_prediction_algorithms(self):
        headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
        potential_algorithms = PredictionClass.prediction_methods()
        prediction_algorithms = []
        for algorithm in potential_algorithms:
            if algorithm in ['NetMHCpanEL', 'NetMHCIIpanEL', 'BigMHC_EL', 'BigMHC_IM', 'DeepImmuno']:
                continue
            if "{} MT IC50 Score".format(algorithm) in headers or "{} IC50 Score".format(algorithm) in headers:
                prediction_algorithms.append(algorithm)
        return prediction_algorithms

    def determine_used_epitope_lengths(self):
        col_name = self.determine_epitope_seq_column_name()
        return list(set([len(s) for s in pd.read_csv(self.input_file, delimiter="\t", usecols=[col_name])[col_name]]))

    def determine_epitope_seq_column_name(self):
        headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
        for header in ["MT Epitope Seq", "Epitope Seq"]:
            if header in headers:
                return header
        raise Exception("No mutant epitope sequence header found.")

    def problematic_positions_exist(self):
        headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
        return 'Problematic Positions' in headers

    def mane_select_run(self):
        df = pd.read_csv(self.input_file, delimiter="\t", nrows=1)
        return df.iloc[0]['MANE Select'] != 'Not Run'

    def canonical_run(self):
        df = pd.read_csv(self.input_file, delimiter="\t", nrows=1)
        return df.iloc[0]['Canonical'] != 'Not Run'

    def calculate_allele_expr(self, line):
        if line['Gene Expression'] == 'NA' or line['Tumor RNA VAF'] == 'NA':
            return 'NA'
        else:
            return round(float(line['Gene Expression']) * float(line['Tumor RNA VAF']), 3)

    def determine_used_el_algorithms(self):
        headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
        potential_algorithms = ["MHCflurryEL Processing", "MHCflurryEL Presentation", "NetMHCpanEL", "NetMHCIIpanEL", "BigMHC_EL", 'BigMHC_IM', 'DeepImmuno']
        prediction_algorithms = []
        for algorithm in potential_algorithms:
            if "{} MT Score".format(algorithm) in headers or "{} Score".format(algorithm) in headers:
                prediction_algorithms.append(algorithm)
        return prediction_algorithms

    def determine_used_percentile_algorithms(self, prediction_algorithms, el_algorithms):
        headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
        percentile_algorithms = []
        for algorithm in prediction_algorithms + el_algorithms:
            if "{} MT Percentile".format(algorithm) in headers or "{} Percentile".format(algorithm) in headers:
                percentile_algorithms.append(algorithm)
        return percentile_algorithms

    def determine_columns_used_for_aggregation(self, prediction_algorithms, el_algorithms):
        used_columns = [
            "Index", "Chromosome", "Start", "Stop", "Reference", "Variant",
            "Transcript", "Transcript Support Level", "MANE Select", "Canonical", "Biotype", "Transcript CDS Flags", "Transcript Length", "Variant Type", "Mutation",
            "Protein Position", "Gene Name", "HLA Allele",
            "Mutation Position", "MT Epitope Seq", "WT Epitope Seq",
            "Tumor DNA VAF", "Tumor RNA Depth",
            "Tumor RNA VAF", "Gene Expression", "Transcript Expression",
            "Median MT IC50 Score", "Median WT IC50 Score", "Median MT Percentile", "Median WT Percentile",
            "Best MT IC50 Score", "Corresponding WT IC50 Score", "Best MT Percentile", "Corresponding WT Percentile",
        ]
        for algorithm in prediction_algorithms:
            used_columns.extend(["{} WT IC50 Score".format(algorithm), "{} MT IC50 Score".format(algorithm)])
            used_columns.extend(["{} WT Percentile".format(algorithm), "{} MT Percentile".format(algorithm)])
        for algorithm in el_algorithms:
            used_columns.extend(["{} WT Score".format(algorithm), "{} MT Score".format(algorithm)])
            if algorithm not in ["MHCflurryEL Processing", "BigMHC_EL", "BigMHC_IM", 'DeepImmuno']:
                used_columns.extend(["{} WT Percentile".format(algorithm), "{} MT Percentile".format(algorithm)])
        if self.problematic_positions_exist():
            used_columns.append("Problematic Positions")
        return used_columns

    def set_column_types(self, prediction_algorithms):
        dtypes = {
            'Chromosome': "string",
            "Start": "int32",
            "Stop": "int32",
            'Reference': "string",
            'Variant': "string",
            "Variant Type": "category",
            "Mutation Position": "category",
            "Median MT IC50 Score": "float32",
            "Median MT Percentile": "float32",
            "Best MT IC50 Score": "float32",
            "Best MT Percentile": "float32",
            "Protein Position": "string",
            "Transcript Length": "int32",
        }
        if self.__class__.__name__ == 'PvacseqAggregateAllEpitopes':
            if self.mane_select_run():
                dtypes['MANE Select'] = 'boolean'
            if self.canonical_run():
                dtypes['Canonical'] = 'boolean'
        for algorithm in prediction_algorithms:
            if algorithm == 'SMM' or algorithm == 'SMMPMBEC':
                continue
            dtypes["{} MT Score".format(algorithm)] = "float32"
            dtypes["{} MT Percentile".format(algorithm)] = "float32"
        return dtypes

    def execute(self):
        prediction_algorithms = self.determine_used_prediction_algorithms()
        epitope_lengths = self.determine_used_epitope_lengths()
        el_algorithms = self.determine_used_el_algorithms()
        percentile_algorithms = self.determine_used_percentile_algorithms(prediction_algorithms, el_algorithms)
        used_columns = self.determine_columns_used_for_aggregation(prediction_algorithms, el_algorithms)
        dtypes = self.set_column_types(prediction_algorithms)

        ##do a crude estimate of clonal vaf/purity
        self.calculate_clonal_vaf()

        if self.__class__.__name__ == 'PvacseqAggregateAllEpitopes':
            metrics = {
                'tumor_purity': self.tumor_purity,
                'vaf_clonal': round(self.vaf_clonal, 3),
                'vaf_subclonal': round(self.vaf_clonal/2, 3),
                'binding_threshold': self.binding_threshold,
                'aggregate_inclusion_binding_threshold': self.aggregate_inclusion_binding_threshold,
                'aggregate_inclusion_count_limit': self.aggregate_inclusion_count_limit,
                'trna_vaf': self.trna_vaf,
                'trna_cov': self.trna_cov,
                'allele_expr_threshold': self.allele_expr_threshold,
                'transcript_prioritization_strategy': self.transcript_prioritization_strategy,
                'maximum_transcript_support_level': self.maximum_transcript_support_level,
                'percentile_threshold': self.percentile_threshold,
                'percentile_threshold_strategy': self.percentile_threshold_strategy,
                'use_allele_specific_binding_thresholds': self.use_allele_specific_binding_thresholds,
                'mt_top_score_metric': self.mt_top_score_metric,
                'wt_top_score_metric': self.wt_top_score_metric,
                'top_score_metric2': self.top_score_metric2,
                'allele_specific_binding_thresholds': self.allele_specific_binding_thresholds,
                'allele_specific_anchors': self.anchor_calculator.use_allele_specific_anchors,
                'alleles': self.hla_types.tolist(),
                'anchor_contribution_threshold': self.anchor_calculator.anchor_contribution_threshold,
                'epitope_lengths': epitope_lengths,
            }
        else:
            metrics = {}

        data = []
        all_epitopes_df = self.read_input_file(used_columns, dtypes)

        ## get a list of unique mutations
        keys = self.get_list_unique_mutation_keys(all_epitopes_df)

        for key in keys:
            (df, key_str) = self.get_sub_df(all_epitopes_df, key)
            (best_mut_line, metrics_for_key) = self.get_best_mut_line(df, key_str, prediction_algorithms, el_algorithms, percentile_algorithms)
            data.append(best_mut_line)
            metrics[key_str] = metrics_for_key
        peptide_table = pd.DataFrame(data=data)

        peptide_table.to_csv(self.output_file, sep='\t', na_rep='NA', index=False, float_format='%.3f')

        self.tier_aggregated_report()

        self.write_metrics_file(metrics)
        self.copy_pvacview_r_files()


class PvacseqAggregateAllEpitopes(AggregateAllEpitopes, metaclass=ABCMeta):
    def __init__(
            self,
            input_file,
            output_file,
            tumor_purity=None,
            binding_threshold=500,
            trna_vaf=0.25,
            trna_cov=10,
            expn_val=1,
            transcript_prioritization_strategy=['canonical', 'mane_select', 'tsl'],
            maximum_transcript_support_level=1,
            percentile_threshold=None,
            percentile_threshold_strategy='conservative',
            allele_specific_binding_thresholds=False,
            top_score_metric="median",
            top_score_metric2="ic50",
            allele_specific_anchors=False,
            allow_incomplete_transcripts=False,
            anchor_contribution_threshold=0.8,
            aggregate_inclusion_binding_threshold=5000,
            aggregate_inclusion_count_limit=15,
        ):
        self.input_file = input_file
        self.output_file = output_file
        self.tumor_purity = tumor_purity
        self.binding_threshold = binding_threshold
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
        self.percentile_threshold = percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.aggregate_inclusion_binding_threshold = aggregate_inclusion_binding_threshold
        self.aggregate_inclusion_count_limit = aggregate_inclusion_count_limit
        self.expn_val = expn_val
        self.allele_expr_threshold = trna_vaf * expn_val * 10
        self.trna_cov = trna_cov
        self.trna_vaf = trna_vaf
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level
        self.allow_incomplete_transcripts = allow_incomplete_transcripts
        if top_score_metric == 'median':
            self.mt_top_score_metric = "Median"
            self.wt_top_score_metric = "Median"
        else:
            self.mt_top_score_metric = "Best"
            self.wt_top_score_metric = "Corresponding"
        self.top_score_metric2 = top_score_metric2
        if self.top_score_metric2 == "percentile":
            self.top_score_mode = "Percentile"
        else:
            self.top_score_mode = "IC50 Score"
        self.metrics_file = output_file.replace('.tsv', '.metrics.json')
        super().__init__()
        self.anchor_calculator = AnchorResiduePass(binding_threshold, self.use_allele_specific_binding_thresholds, self.allele_specific_binding_thresholds, allele_specific_anchors, anchor_contribution_threshold, self.wt_top_score_metric)

    def get_list_unique_mutation_keys(self, df):
        keys = df[['Chromosome', 'Start', 'Stop', 'Reference', 'Variant']].values.tolist()
        keys = [list(i) for i in set(tuple(i) for i in keys)]
        return sorted(keys)

    def calculate_clonal_vaf(self):
        if self.tumor_purity:
            self.vaf_clonal =  self.tumor_purity * 0.5
            print("Tumor clonal VAF estimated as {} (calculated from user-provided tumor purity of {}). Assuming variants with VAF < {} are subclonal".format(round(self.vaf_clonal, 3), round(self.tumor_purity, 3), round(self.vaf_clonal/2, 3)))
        else:
        #if no tumor purity is provided, make a rough estimate by taking the list of VAFs < 0.6 (assumption is that these are CN-neutral) and return the largest as the marker of the founding clone
            vafs = np.sort(pd.read_csv(self.input_file, delimiter="\t", usecols=["Tumor DNA VAF"])['Tumor DNA VAF'].unique())[::-1]
            vafs_clonal = list(filter(lambda vaf: vaf < 0.6, vafs))
            if len(vafs_clonal) == 0:
                self.vaf_clonal = 0.6
            else:
                self.vaf_clonal = vafs_clonal[0]
                if self.vaf_clonal > 0.5:
                    self.vaf_clonal = 0.5
            print("Tumor clonal VAF estimated as {} (estimated from Tumor DNA VAF data). Assuming variants with VAF < {} are subclonal".format(round(self.vaf_clonal, 3), round(self.vaf_clonal/2, 3)))
            return self.vaf_clonal

    def read_input_file(self, used_columns, dtypes):
        df = pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False, na_values="NA", keep_default_na=False, usecols=used_columns, dtype=dtypes)
        df = df.astype({"{} MT IC50 Score".format(self.mt_top_score_metric):'float'})
        return df

    def get_sub_df(self, all_epitopes_df, key):
        key_str = "{}-{}-{}-{}-{}".format(key[0], key[1], key[2], key[3], key[4])
        df = (all_epitopes_df[lambda x: (x['Chromosome'] == key[0]) & (x['Start'] == key[1]) & (x['Stop'] == key[2]) & (x['Reference'] == key[3]) & (x['Variant'] == key[4])]).copy()
        df['annotation'] = df[['Transcript', 'Gene Name', 'Mutation', 'Protein Position']].agg('-'.join, axis=1)
        df['key'] = key_str
        return (df, key_str)

    def get_best_binder(self, df):
        return PvacseqBestCandidate(
            self.transcript_prioritization_strategy,
            self.maximum_transcript_support_level,
            self.anchor_calculator,
            self.mt_top_score_metric,
            self.top_score_mode,
            self.allow_incomplete_transcripts,
        ).get(df)

    def get_included_df(self, df):
        binding_df = df[df["{} MT IC50 Score".format(self.mt_top_score_metric)] < self.aggregate_inclusion_binding_threshold]
        if binding_df.shape[0] == 0:
            return binding_df
        else:
            peptides = list(set(binding_df["MT Epitope Seq"]))
            if len(peptides) <= self.aggregate_inclusion_count_limit:
                return binding_df

            best_peptide_entries = []
            for peptide in peptides:
                peptide_df = binding_df[binding_df["MT Epitope Seq"] == peptide]
                best_peptide_entries.append(self.get_best_binder(peptide_df))
            best_peptide_entries_df = pd.DataFrame(best_peptide_entries)
            top_n_best_peptide_entries_df = self.sort_included_df(best_peptide_entries_df).iloc[:self.aggregate_inclusion_count_limit]
            top_n_best_peptides = list(set(top_n_best_peptide_entries_df["MT Epitope Seq"]))
            return binding_df[binding_df["MT Epitope Seq"].isin(top_n_best_peptides)]

    def get_good_binders_df(self, included_df):
        if self.use_allele_specific_binding_thresholds:
            selection = []
            for index, row in included_df.iterrows():
                if row['HLA Allele'] in self.allele_specific_binding_thresholds:
                    binding_threshold = self.allele_specific_binding_thresholds[row['HLA Allele']]
                else:
                    binding_threshold = self.binding_threshold
                if row["{} MT IC50 Score".format(self.mt_top_score_metric)] < binding_threshold:
                    selection.append(index)
            good_binders = included_df[included_df.index.isin(selection)]
        else:
            good_binders = included_df[included_df["{} MT IC50 Score".format(self.mt_top_score_metric)] < self.binding_threshold]
        return good_binders

    def sort_included_df(self, df):
        df['biotype_sort'] = df['Biotype'].apply(lambda x: 1 if x == 'protein_coding' else 2)
        df['tsl_sort'] = df['Transcript Support Level'].apply(lambda x: 6 if pd.isnull(x) or x == 'Not Supported' else int(x))
        df['anchor_residue_pass_sort'] = df.apply(lambda x: 1 if self.anchor_calculator.is_anchor_residue_pass(x) else 2, axis=1)

        sort_columns = ["biotype_sort", "tsl_sort", "anchor_residue_pass_sort"]
        sort_order = [True, True, True]

        if self.allow_incomplete_transcripts:
            df['transcript_cds_flags_sort'] = df['Transcript CDS Flags'].apply(lambda x: 1 if x == "None" else (2 if any(flag in str(x) for flag in ["cds_start_nf", "cds_end_nf"]) else 1))
            sort_columns.append('transcript_cds_flags_sort')
            sort_order.append(True)

        if self.problematic_positions_exist():
            df['problematic_positions_sort'] = df['Problematic Positions'].apply(lambda x: 1 if x == "None" else 2)
            sort_columns.append('problematic_positions_sort')
            sort_order.append(True)

        if self.top_score_mode == 'IC50 Score':
            primary = "{} MT IC50 Score".format(self.mt_top_score_metric)
            secondary = "{} MT Percentile".format(self.mt_top_score_metric)
        else:
            primary = "{} MT Percentile".format(self.mt_top_score_metric)
            secondary = "{} MT IC50 Score".format(self.mt_top_score_metric)
        sort_columns.extend([primary, "Transcript Length", secondary])
        sort_order.extend([True, False, True])
        df.sort_values(by=sort_columns, inplace=True, ascending=sort_order)
        df.drop(columns=[
            'biotype_sort',
            'tsl_sort',
            'anchor_residue_pass_sort',
            'transcript_cds_flags_sort',
            'problematic_positions_sort'
        ], inplace=True, errors='ignore')
        return df

    def get_unique_peptide_hla_counts(self, good_binders_df):
        return pd.DataFrame(good_binders_df.groupby(['HLA Allele', 'MT Epitope Seq']).size().reset_index())

    def replace_nas(self, items):
        return ["NA" if pd.isna(x) else x for x in items]

    def round_to_ints(self, items):
        return [round(x) if (type(x) == float and not pd.isna(x)) else x for x in items]

    def get_included_df_metrics(self, included_df, prediction_algorithms, el_algorithms, percentile_algorithms):
        peptides = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        included_peptides = included_df["MT Epitope Seq"].unique()
        included_transcripts = included_df['annotation'].unique()
        peptide_sets = {}
        for annotation in included_transcripts:
            included_df_annotation = included_df[included_df['annotation'] == annotation]
            peptide_set = tuple(included_df_annotation["MT Epitope Seq"].unique())
            if peptide_set in peptide_sets:
                peptide_sets[peptide_set].append(annotation)
            else:
                peptide_sets[peptide_set] = [annotation]

        set_number = 1
        for peptide_set, annotations in peptide_sets.items():
            set_name = "Transcript Set {}".format(set_number)
            annotation = annotations[0]
            included_df_annotation = included_df[included_df['annotation'] == annotation]
            results = defaultdict(lambda: defaultdict(list))
            for peptide in list(peptide_set):
                included_df_peptide_annotation = included_df_annotation[included_df_annotation['MT Epitope Seq'] == peptide]
                if len(included_df_peptide_annotation) > 0:
                    individual_ic50_calls = { 'algorithms': prediction_algorithms }
                    individual_ic50_percentile_calls = { 'algorithms': prediction_algorithms }
                    individual_el_calls = { 'algorithms': el_algorithms }
                    individual_el_percentile_calls = { 'algorithms': el_algorithms }
                    individual_percentile_calls = { 'algorithms': percentile_algorithms }
                    anchor_fails = []
                    for peptide_type, top_score_metric in zip(['MT', 'WT'], [self.mt_top_score_metric, self.wt_top_score_metric]):
                        ic50s = {}
                        percentiles = {}
                        ic50_calls = {}
                        percentile_calls = {}
                        el_calls = {}
                        el_percentile_calls = {}
                        all_percentile_calls = {}
                        for index, line in included_df_peptide_annotation.to_dict(orient='index').items():
                            ic50s[line['HLA Allele']] = line['{} {} IC50 Score'.format(top_score_metric, peptide_type)]
                            percentiles[line['HLA Allele']] = line['{} {} Percentile'.format(top_score_metric, peptide_type)]
                            ic50_calls[line['HLA Allele']] = self.replace_nas([line["{} {} IC50 Score".format(algorithm, peptide_type)] for algorithm in prediction_algorithms])
                            percentile_calls[line['HLA Allele']] = self.replace_nas([line["{} {} Percentile".format(algorithm, peptide_type)] for algorithm in prediction_algorithms])
                            el_calls[line['HLA Allele']] = self.replace_nas([line["{} {} Score".format(algorithm, peptide_type)] for algorithm in el_algorithms])
                            el_percentile_calls[line['HLA Allele']] = self.replace_nas(['NA' if algorithm in ['MHCflurryEL Processing', 'BigMHC_EL', 'BigMHC_IM', 'DeepImmuno'] else line["{} {} Percentile".format(algorithm, peptide_type)] for algorithm in el_algorithms])
                            all_percentile_calls[line['HLA Allele']] = self.replace_nas([line["{} {} Percentile".format(algorithm, peptide_type)] for algorithm in percentile_algorithms])
                            if peptide_type == 'MT' and not self.anchor_calculator.is_anchor_residue_pass(line):
                                anchor_fails.append(line['HLA Allele'])
                        sorted_ic50s = []
                        sorted_percentiles = []
                        for hla_type in sorted(self.hla_types):
                            if hla_type in ic50s:
                                sorted_ic50s.append(ic50s[hla_type])
                            else:
                                sorted_ic50s.append('X')
                            if hla_type in percentiles:
                                sorted_percentiles.append(percentiles[hla_type])
                            else:
                                sorted_percentiles.append('X')
                        results[peptide]['ic50s_{}'.format(peptide_type)] = self.replace_nas(sorted_ic50s)
                        results[peptide]['percentiles_{}'.format(peptide_type)] = self.replace_nas(sorted_percentiles)
                        individual_ic50_calls[peptide_type] = ic50_calls
                        individual_ic50_percentile_calls[peptide_type] = percentile_calls
                        individual_el_calls[peptide_type] = el_calls
                        individual_el_percentile_calls[peptide_type] = el_percentile_calls
                        individual_percentile_calls[peptide_type] = all_percentile_calls
                    results[peptide]['hla_types'] = sorted(self.hla_types)
                    results[peptide]['mutation_position'] = "NA" if pd.isna(included_df_peptide_annotation.iloc[0]['Mutation Position']) else str(included_df_peptide_annotation.iloc[0]['Mutation Position'])
                    results[peptide]['problematic_positions'] = str(included_df_peptide_annotation.iloc[0]['Problematic Positions']) if 'Problematic Positions' in included_df_peptide_annotation.iloc[0] else 'None'
                    if len(anchor_fails) > 0:
                        results[peptide]['anchor_fails'] = ', '.join(anchor_fails)
                    else:
                        results[peptide]['anchor_fails'] = 'None'
                    results[peptide]['individual_ic50_calls'] = individual_ic50_calls
                    results[peptide]['individual_ic50_percentile_calls'] = individual_ic50_percentile_calls
                    results[peptide]['individual_el_calls'] = individual_el_calls
                    results[peptide]['individual_el_percentile_calls'] = individual_el_percentile_calls
                    results[peptide]['individual_percentile_calls'] = individual_percentile_calls
                    wt_peptide = included_df_peptide_annotation.iloc[0]['WT Epitope Seq']
                    if pd.isna(wt_peptide):
                        variant_type = included_df_peptide_annotation.iloc[0]['Variant Type']
                        if variant_type == 'FS':
                            wt_peptide = 'FS-NA'
                        elif variant_type == 'inframe_ins':
                            wt_peptide = 'INS-NA'
                        elif variant_type == 'inframe_del':
                            wt_peptide = 'DEL-NA'
                    results[peptide]['wt_peptide'] = wt_peptide
            peptides[set_name]['peptides'] = self.sort_peptides(results)
            sorted_transcripts = self.sort_transcripts(annotations, included_df)
            peptides[set_name]['transcripts'] = list(sorted_transcripts.Annotation)
            peptides[set_name]['transcript_expr'] = self.replace_nas(list(sorted_transcripts.Expr))
            peptides[set_name]['mane_select'] = list(sorted_transcripts['MANE Select'])
            peptides[set_name]['canonical'] = list(sorted_transcripts.Canonical)
            peptides[set_name]['tsl'] = self.replace_nas(self.round_to_ints(list(sorted_transcripts['Transcript Support Level'])))
            peptides[set_name]['biotype'] = list(sorted_transcripts.Biotype)
            peptides[set_name]['transcript_cds_flags'] = list(sorted_transcripts['Transcript CDS Flags'])
            peptides[set_name]['transcript_length'] = [int(l) for l in list(sorted_transcripts.Length)]
            peptides[set_name]['transcript_count'] = len(annotations)
            peptides[set_name]['peptide_count'] = len(peptide_set)
            peptides[set_name]['total_expr'] = math.fsum([0 if x == 'NA' else (float(x)) for x in peptides[set_name]['transcript_expr']])
            set_number += 1
        anno_count = len(included_transcripts)

        return (peptides, anno_count)

    def sort_peptides(self, results):
        for k, v in results.items():
            v['problematic_positions_sort'] = 1 if v['problematic_positions'] == 'None' else 2
            v['anchor_fail_sort'] = 1 if v['anchor_fails'] == 'None' else 2
            v['best_ic50s_MT'] = min([ic50 for ic50 in v['ic50s_MT'] if ic50 != 'X'])
        sorted_results = dict(sorted(results.items(), key=lambda x:(x[1]['problematic_positions_sort'],x[1]['anchor_fail_sort'],x[1]['best_ic50s_MT'])))
        for k, v in sorted_results.items():
            v.pop('problematic_positions_sort')
            v.pop('anchor_fail_sort')
            v.pop('best_ic50s_MT')
        return sorted_results

    def sort_transcripts(self, annotations, included_df):
        records = []
        for annotation in annotations:
            line = included_df[included_df['annotation'] == annotation].iloc[0]
            records.append({
                'Annotation': line['annotation'],
                'Biotype': line['Biotype'],
                'Transcript CDS Flags': line['Transcript CDS Flags'],
                'MANE Select': line['MANE Select'],
                'Canonical': line['Canonical'],
                'Transcript Support Level': line['Transcript Support Level'],
                'Length': line['Transcript Length'],
                'Expr': line['Transcript Expression'],
            })
        transcript_table = pd.DataFrame.from_records(records)

        transcript_table['Biotype Sort'] = transcript_table.Biotype.map(lambda x: 1 if x == 'protein_coding' else 2)
        transcript_table['TSL Sort'] = transcript_table.apply(lambda x: 1 if is_preferred_transcript(x, self.transcript_prioritization_strategy, self.maximum_transcript_support_level) else 2, axis=1)

        if self.allow_incomplete_transcripts:
            transcript_table['Transcript CDS Flags Sort'] = transcript_table['Transcript CDS Flags'].apply(lambda x: 1 if x == "None" else (2 if any(flag in str(x) for flag in ["cds_start_nf", "cds_end_nf"]) else 1))
            sort_columns = ["Biotype Sort", "Transcript CDS Flags Sort", "TSL Sort", "Length"]
            sort_order = [True, True, True, False]
        else:
            sort_columns = ["Biotype Sort", "TSL Sort", "Length"]
            sort_order = [True, True, False]

        transcript_table.sort_values(by=sort_columns, ascending=sort_order, inplace=True)
        return transcript_table

    def calculate_unique_peptide_count(self, included_df):
        return len(included_df["MT Epitope Seq"].unique())

    def calculate_good_binder_count(self, good_binders_df):
        return len(good_binders_df["MT Epitope Seq"].unique())

    def get_default_annotation_count(self):
        return 0

    def get_best_aa_change(self, best):
        if best['Variant Type'] == 'FS':
            return 'FS{}'.format(best['Protein Position'])
        else:
            (wt_aa, mt_aa) = best["Mutation"].split("/")
            return "".join([wt_aa, best["Protein Position"], mt_aa])

    def assemble_result_line(self, best, key, hla, anno_count, included_peptide_count, good_binder_count):
        allele_expr = self.calculate_allele_expr(best)

        problematic_positions = best['Problematic Positions'] if 'Problematic Positions' in best else 'None'

        tsl = best['Transcript Support Level'] if best['Transcript Support Level'] == "Not Supported" or pd.isna(best['Transcript Support Level']) else str(int(best['Transcript Support Level']))

        out_dict = { 'ID': key, 'Index': best['Index'] }
        out_dict.update({ k.replace('HLA-', ''):v for k,v in sorted(hla.items()) })
        out_dict.update({
            'Gene': best["Gene Name"],
            'AA Change': self.get_best_aa_change(best),
            'Num Passing Transcripts': anno_count,
            'Best Peptide': best["MT Epitope Seq"],
            'Best Transcript': best["Transcript"],
            'MANE Select': best['MANE Select'],
            'Canonical': best['Canonical'],
            'TSL': tsl,
            'Allele': best["HLA Allele"],
            'Pos': best["Mutation Position"],
            'Prob Pos': problematic_positions,
            'Num Included Peptides': included_peptide_count,
            'Num Passing Peptides': good_binder_count,
            'IC50 MT': best["{} MT IC50 Score".format(self.mt_top_score_metric)],
            'IC50 WT': best["{} WT IC50 Score".format(self.wt_top_score_metric)],
            '%ile MT': best["{} MT Percentile".format(self.mt_top_score_metric)],
            '%ile WT': best["{} WT Percentile".format(self.wt_top_score_metric)],
            'RNA Expr': best["Gene Expression"],
            'RNA VAF': best["Tumor RNA VAF"],
            'Allele Expr': allele_expr,
            'RNA Depth': best["Tumor RNA Depth"],
            'DNA VAF': best["Tumor DNA VAF"],
            'Tier': 'Not Tiered',
            'Evaluation': 'Pending',
        })
        return out_dict

    def get_metrics(self, peptides, best):
        return {
            'good_binders': peptides,
            'sets': list(peptides.keys()),
            'transcript_counts': [v['transcript_count'] for k, v in peptides.items()],
            'peptide_counts': [v['peptide_count'] for k, v in peptides.items()],
            'set_expr': [v['total_expr'] for k, v in peptides.items()],
            'DNA VAF': 'NA' if pd.isna(best['Tumor DNA VAF']) else float(best['Tumor DNA VAF']),
            'RNA VAF': 'NA' if pd.isna(best['Tumor RNA VAF']) else float(best['Tumor RNA VAF']),
            'gene_expr': 'NA' if pd.isna(best['Gene Expression']) else float(best['Gene Expression']),
            'best_peptide_mt': best['MT Epitope Seq'],
            'best_peptide_wt': 'NA' if pd.isna(best['WT Epitope Seq']) else best['WT Epitope Seq'],
            'best_hla_allele': best['HLA Allele'],
        }

    def write_metrics_file(self, metrics):
        with open(self.metrics_file, 'w') as fh:
            json.dump(metrics, fh, indent=2, separators=(',', ': '))

    def copy_pvacview_r_files(self):
        module_dir = os.path.dirname(__file__)
        r_folder = os.path.abspath(os.path.join(module_dir,"..","tools","pvacview"))
        files = glob.iglob(os.path.join(r_folder, "*.R"))
        destination = os.path.abspath(os.path.dirname(self.output_file))
        os.makedirs(os.path.join(destination, "www"), exist_ok=True)
        for i in files:
            shutil.copy(i, destination)
        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            shutil.copy(os.path.join(r_folder, "www", i), os.path.join(destination, "www", i))

    def tier_aggregated_report(self):
        PvacseqUpdateTiers(
            self.output_file,
            self.vaf_clonal,
            binding_threshold=self.binding_threshold,
            trna_vaf=self.trna_vaf,
            trna_cov=self.trna_cov,
            expn_val=self.expn_val,
            transcript_prioritization_strategy=self.transcript_prioritization_strategy,
            maximum_transcript_support_level=self.maximum_transcript_support_level,
            percentile_threshold=self.percentile_threshold,
            percentile_threshold_strategy=self.percentile_threshold_strategy,
            allele_specific_binding_thresholds=self.use_allele_specific_binding_thresholds,
            allele_specific_anchors=self.anchor_calculator.use_allele_specific_anchors,
            anchor_contribution_threshold=self.anchor_calculator.anchor_contribution_threshold,
            top_score_metric2=self.top_score_metric2,
        ).execute()


class UnmatchedSequenceAggregateAllEpitopes(AggregateAllEpitopes, metaclass=ABCMeta):
    def __init__(self,
            input_file,
            output_file,
            binding_threshold=500,
            percentile_threshold=None,
            percentile_threshold_strategy='conservative',
            allele_specific_binding_thresholds=False,
            top_score_metric="median",
            top_score_metric2="ic50",
            aggregate_inclusion_binding_threshold=5000,
            aggregate_inclusion_count_limit=15,
        ):
        self.input_file = input_file
        self.output_file = output_file
        self.binding_threshold = binding_threshold
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
        self.percentile_threshold = percentile_threshold
        self.percentile_threshold_strategy = percentile_threshold_strategy
        self.aggregate_inclusion_binding_threshold = aggregate_inclusion_binding_threshold
        self.aggregate_inclusion_count_limit = aggregate_inclusion_count_limit
        if top_score_metric == 'median':
            self.top_score_metric = "Median"
        else:
            self.top_score_metric = "Best"
        self.top_score_metric2 = top_score_metric2
        if self.top_score_metric2 == "percentile":
            self.top_score_mode = "Percentile"
        else:
            self.top_score_mode = "IC50 Score"
        self.metrics_file = output_file.replace('.tsv', '.metrics.json')
        super().__init__()

    def get_list_unique_mutation_keys(self, df):
        keys = df["Mutation"].values.tolist()
        return sorted(list(set(keys)))

    def calculate_clonal_vaf(self):
        if self.__class__.__name__ == 'PvacspliceAggregateAllEpitopes':
            PvacseqAggregateAllEpitopes.calculate_clonal_vaf(self)
        else:
            self.vaf_clonal = None

    def read_input_file(self, used_columns, dtypes):
        df = pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False, na_values="NA", keep_default_na=False, dtype={"Mutation": str})
        df = df[df["{} IC50 Score".format(self.top_score_metric)] != 'NA']
        df = df.astype({"{} IC50 Score".format(self.top_score_metric):'float'})
        return df

    def get_sub_df(self, all_epitopes_df, key):
        df = (all_epitopes_df[lambda x: (x['Mutation'] == key)]).copy()
        return (df, key)

    def get_included_df(self, df):
        binding_df = df[df[f"{self.top_score_metric} IC50 Score"] < self.aggregate_inclusion_binding_threshold]

        if binding_df.shape[0] == 0:
            return binding_df
        else:
            peptides = list(set(binding_df["Epitope Seq"]))
            if len(peptides) <= self.aggregate_inclusion_count_limit:
                return binding_df

            best_peptide_entries = []
            for peptide in peptides:
                peptide_df = binding_df[binding_df["Epitope Seq"] == peptide]
                best_peptide_entries.append(self.get_best_binder(peptide_df))
            best_peptide_entries_df = pd.DataFrame(best_peptide_entries)
            top_n_best_peptide_entries_df = self.sort_included_df(best_peptide_entries_df).iloc[:self.aggregate_inclusion_count_limit]
            top_n_best_peptides = list(set(top_n_best_peptide_entries_df["Epitope Seq"]))
            return binding_df[binding_df["Epitope Seq"].isin(top_n_best_peptides)]

    def get_good_binders_df(self, included_df):
        if self.use_allele_specific_binding_thresholds:
            selection = []
            for index, row in included_df.iterrows():
                if row['HLA Allele'] in self.allele_specific_binding_thresholds:
                    binding_threshold = self.allele_specific_binding_thresholds[row['HLA Allele']]
                else:
                    binding_threshold = self.binding_threshold
                if row["{} IC50 Score".format(self.top_score_metric)] < binding_threshold:
                    selection.append(index)
            good_binders = included_df[included_df.index.isin(selection)]
        else:
            good_binders = included_df[included_df["{} IC50 Score".format(self.top_score_metric)] < self.binding_threshold]
        return good_binders

    def sort_included_df(self, df):
        if self.top_score_metric2 == "percentile":
            df.sort_values(by=[
                "{} Percentile".format(self.top_score_metric),
                "{} IC50 Score".format(self.top_score_metric),
            ], inplace=True, ascending=[True, True])
        else:
            df.sort_values(by=[
                "{} IC50 Score".format(self.top_score_metric),
                "{} Percentile".format(self.top_score_metric),
            ], inplace=True, ascending=[True, True])
        return df

    def get_unique_peptide_hla_counts(self, good_binders_df):
        return pd.DataFrame(good_binders_df.groupby(['HLA Allele', 'Epitope Seq']).size().reset_index())

    def get_included_df_metrics(self, included_df, prediction_algorithms, el_algorithms, percentile_algorithms):
        return (None, "NA")

    def calculate_unique_peptide_count(self, included_df):
        return len(included_df["Epitope Seq"].unique())

    def calculate_good_binder_count(self, good_binders_df):
        return len(good_binders_df["Epitope Seq"].unique())

    def get_default_annotation_count(self):
        return "NA"

    def get_metrics(self, peptides, best):
        return None

    def write_metrics_file(self, metrics):
        pass

    def copy_pvacview_r_files(self):
        pass

    def tier_aggregated_report(self):
        pass

class PvacfuseAggregateAllEpitopes(UnmatchedSequenceAggregateAllEpitopes, metaclass=ABCMeta):
    def __init__(
        self,
        input_file,
        output_file,
        binding_threshold=500,
        percentile_threshold=None,
        percentile_threshold_strategy='conservative',
        allele_specific_binding_thresholds=False,
        top_score_metric="median",
        top_score_metric2="ic50",
        read_support=5,
        expn_val=0.1,
        aggregate_inclusion_binding_threshold=5000,
        aggregate_inclusion_count_limit=15,
    ):
        UnmatchedSequenceAggregateAllEpitopes.__init__(
            self,
            input_file,
            output_file,
            binding_threshold=binding_threshold,
            percentile_threshold=percentile_threshold,
            percentile_threshold_strategy = percentile_threshold_strategy,
            allele_specific_binding_thresholds=allele_specific_binding_thresholds,
            top_score_metric=top_score_metric,
            top_score_metric2=top_score_metric2,
            aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
            aggregate_inclusion_count_limit=aggregate_inclusion_count_limit,
        )
        self.read_support = read_support
        self.expn_val = expn_val

    def assemble_result_line(self, best, key, hla, anno_count, included_peptide_count, good_binder_count):
        out_dict = { 'ID': key }
        out_dict.update({ k.replace('HLA-', ''):v for k,v in sorted(hla.items()) })
        gene = best['Gene Name'] if 'Gene Name' in best else 'NA'
        transcript = best['Transcript'] if 'Transcript' in best else 'NA'
        problematic_positions = best['Problematic Positions'] if 'Problematic Positions' in best else 'None'
        out_dict.update({
            'Gene': gene,
            'Best Peptide': best["Epitope Seq"],
            'Best Transcript': transcript,
            'Allele': best['HLA Allele'],
            'Prob Pos': problematic_positions,
            'Num Included Peptides': included_peptide_count,
            'Num Passing Peptides': good_binder_count,
            'IC50 MT': best["{} IC50 Score".format(self.top_score_metric)],
            '%ile MT': best["{} Percentile".format(self.top_score_metric)],
            'Expr': best['Expression'],
            'Read Support': best['Read Support'],
            'Tier': 'Not Tiered',
            'Evaluation': 'Pending',
        })
        return out_dict

    def get_best_binder(self, df):
        return PvacfuseBestCandidate(
            self.top_score_metric,
            self.top_score_mode,
        ).get(df)

    def tier_aggregated_report(self):
        PvacfuseUpdateTiers(
            self.output_file,
            binding_threshold=self.binding_threshold,
            allele_specific_binding_thresholds=self.use_allele_specific_binding_thresholds,
            percentile_threshold=self.percentile_threshold,
            percentile_threshold_strategy=self.percentile_threshold_strategy,
            read_support=self.read_support,
            expn_val=self.expn_val,
            top_score_metric2=self.top_score_metric2,
        ).execute()


class PvacbindAggregateAllEpitopes(UnmatchedSequenceAggregateAllEpitopes, metaclass=ABCMeta):
    def assemble_result_line(self, best, key, hla, anno_count, included_peptide_count, good_binder_count):
        out_dict = { 'ID': key }
        out_dict.update({ k.replace('HLA-', ''):v for k,v in sorted(hla.items()) })
        problematic_positions = best['Problematic Positions'] if 'Problematic Positions' in best else 'None'
        out_dict.update({
            'Best Peptide': best["Epitope Seq"],
            'Allele': best["HLA Allele"],
            'Pos': best['Sub-peptide Position'],
            'Prob Pos': problematic_positions,
            'Num Included Peptides': included_peptide_count,
            'Num Passing Peptides': good_binder_count,
            'IC50 MT': best["{} IC50 Score".format(self.top_score_metric)],
            '%ile MT': best["{} Percentile".format(self.top_score_metric)],
            'Tier': 'Not Tiered',
            'Evaluation': 'Pending',
        })
        return out_dict

    def get_best_binder(self, df):
        return PvacbindBestCandidate(
            self.top_score_metric,
            self.top_score_mode,
        ).get(df)

    def tier_aggregated_report(self):
        PvacbindUpdateTiers(
            self.output_file,
            binding_threshold=self.binding_threshold,
            allele_specific_binding_thresholds=self.use_allele_specific_binding_thresholds,
            percentile_threshold=self.percentile_threshold,
            percentile_threshold_strategy=self.percentile_threshold_strategy,
            top_score_metric2=self.top_score_metric2,
        ).execute()


class PvacspliceAggregateAllEpitopes(PvacbindAggregateAllEpitopes, metaclass=ABCMeta):
    def __init__(
        self,
        input_file,
        output_file,
        tumor_purity=None,
        binding_threshold=500,
        percentile_threshold=None,
        percentile_threshold_strategy='conservative',
        allele_specific_binding_thresholds=False,
        aggregate_inclusion_binding_threshold=5000,
        aggregate_inclusion_count_limit=15,
        top_score_metric="median",
        top_score_metric2="ic50",
        trna_vaf=0.25,
        trna_cov=10,
        expn_val=1,
        transcript_prioritization_strategy=['canonical', 'mane_select', 'tsl'],
        maximum_transcript_support_level=1,
        allow_incomplete_transcripts=False,
    ):
        PvacbindAggregateAllEpitopes.__init__(
            self,
            input_file,
            output_file,
            binding_threshold=binding_threshold,
            percentile_threshold=percentile_threshold,
            percentile_threshold_strategy = percentile_threshold_strategy,
            allele_specific_binding_thresholds=allele_specific_binding_thresholds,
            aggregate_inclusion_binding_threshold=aggregate_inclusion_binding_threshold,
            aggregate_inclusion_count_limit=aggregate_inclusion_count_limit,
            top_score_metric=top_score_metric,
            top_score_metric2=top_score_metric2,
        )
        self.tumor_purity = tumor_purity
        self.trna_vaf = trna_vaf
        self.trna_cov = trna_cov
        self.expn_val = expn_val
        self.allele_expr_threshold = trna_vaf * expn_val * 10
        self.transcript_prioritization_strategy = transcript_prioritization_strategy
        self.maximum_transcript_support_level = maximum_transcript_support_level
        self.allow_incomplete_transcripts = allow_incomplete_transcripts

    # pvacbind w/ Index instead of Mutation
    def get_list_unique_mutation_keys(self, df):
        keys = df["Junction"].values.tolist()
        return sorted(list(set(keys)))

    # pvacbind w/ Index instead of Mutation
    def read_input_file(self, used_columns, dtypes):
        return pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False,
                           na_values="NA", keep_default_na=False, dtype={"Index": str})

    # pvacbind w/ Index instead of Mutation
    def get_sub_df(self, all_epitopes_df, df_key):
        df = (all_epitopes_df[lambda x: (x['Junction'] == df_key)]).copy()
        return df, df_key

    # pvacbind w/ vaf and expression info included
    def assemble_result_line(self, best, key, hla, anno_count, included_peptide_count, good_binder_count):
        out_dict = {'ID': best['Index']}
        out_dict.update({k.replace('HLA-', ''): v for k, v in sorted(hla.items())})

        gene = best['Gene Name'] if 'Gene Name' in best else 'NA'
        transcript = best['Transcript'] if 'Transcript' in best else 'NA'
        problematic_positions = best['Problematic Positions'] if 'Problematic Positions' in best else 'None'
        tsl = best['Transcript Support Level'] if best['Transcript Support Level'] == "Not Supported" or pd.isna(best['Transcript Support Level']) else str(int(best['Transcript Support Level']))
        allele_expr = self.calculate_allele_expr(best)

        out_dict.update({
            'Gene': gene,
            'Transcript': transcript,
            'Junction Name': best['Junction'],
            'AA Change': best['Amino Acid Change'],
            'Best Peptide': best["Epitope Seq"],
            'TSL': tsl,
            'MANE Select': best['MANE Select'],
            'Canonical': best['Canonical'],
            'Allele': best["HLA Allele"],
            'Pos': best['Protein Position'],
            'Prob Pos': problematic_positions,
            'Num Included Peptides': included_peptide_count,
            'Num Passing Peptides': good_binder_count,
            'IC50 MT': best["{} IC50 Score".format(self.top_score_metric)],
            '%ile MT': best["{} Percentile".format(self.top_score_metric)],
            'RNA Expr': best["Gene Expression"],
            'RNA VAF': best["Tumor RNA VAF"],
            'Allele Expr': allele_expr,
            'RNA Depth': best["Tumor RNA Depth"],
            'DNA VAF': best["Tumor DNA VAF"],
            'Tier': 'Not Tiered',
            'Evaluation': 'Pending',
        })
        return out_dict

    def get_best_binder(self, df):
        return PvacspliceBestCandidate(
            self.transcript_prioritization_strategy,
            self.maximum_transcript_support_level,
            self.top_score_metric,
            self.top_score_mode,
            self.allow_incomplete_transcripts,
        ).get(df)

    def tier_aggregated_report(self):
        PvacspliceUpdateTiers(
            self.output_file,
            self.vaf_clonal,
            binding_threshold=self.binding_threshold,
            allele_specific_binding_thresholds=self.use_allele_specific_binding_thresholds,
            percentile_threshold=self.percentile_threshold,
            percentile_threshold_strategy=self.percentile_threshold_strategy,
            trna_vaf=self.trna_vaf,
            trna_cov=self.trna_cov,
            expn_val=self.expn_val,
            transcript_prioritization_strategy=self.transcript_prioritization_strategy,
            maximum_transcript_support_level=self.maximum_transcript_support_level,
            top_score_metric2=self.top_score_metric2,
        ).execute()
