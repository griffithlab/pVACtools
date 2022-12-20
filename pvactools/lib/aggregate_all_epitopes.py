import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import json
from Bio import SeqIO
import os
import shutil
from abc import ABCMeta, abstractmethod

from pvactools.lib.prediction_class import PredictionClass

class AggregateAllEpitopes:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
        self.metrics_file = output_file.replace('.tsv', '.metrics.json')


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
    def get_tier(self, mutation, vaf_clonal):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_good_binders(self, df, max_ic50):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_unique_good_binders(self, good_binders):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_good_binders_metrics(self, good_binders, prediction_algorithms, hla_types):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def calculate_annotation_count(self, peptides):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def calculate_unique_peptide_count(self, good_binders):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_default_annotation_count(self):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_best_aa_change(self, best):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def calculate_allele_expr(self, line):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def assemble_result_line(self, best, key, vaf_clonal, hla, anno_count, peptide_count):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def get_metrics(self, df, peptides, best):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def write_metrics_file(self, metrics):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def sort_table(self, df):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def copy_pvacview_r_files(self):
        raise Exception("Must implement method in child class")


    def get_best_mut_line(self, df, key, hla_types, prediction_algorithms, vaf_clonal, max_ic50=1000):
        #order by best median score and get best ic50 peptide
        best = self.get_best_binder(df)

        #these counts should represent only the "good binders" with ic50 < max
        #for all sites other than tier4 slop
        good_binders = self.get_good_binders(df, max_ic50)
        if len(good_binders) > 0:
            good_binders_uniq = self.get_unique_good_binders(good_binders)
            good_binders_hla = Counter(good_binders_uniq["HLA Allele"])
            hla = dict(map(lambda x : (x, good_binders_hla[x]) if x in good_binders_hla else (x, ""), hla_types))
            #get a list of all unique gene/transcript/aa_change combinations
            #store a count of all unique peptides that passed
            peptides = self.get_good_binders_metrics(good_binders, prediction_algorithms, hla_types)
            anno_count = self.calculate_annotation_count(peptides)
            peptide_count = self.calculate_unique_peptide_count(good_binders)
        else:
            hla = dict(map(lambda x : (x, ""), hla_types))
            peptides = {}
            anno_count = self.get_default_annotation_count()
            peptide_count = 0

        #assemble the line
        out_dict = self.assemble_result_line(best, key, vaf_clonal, hla, anno_count, peptide_count);

        metric = self.get_metrics(df, peptides, best)
        return (out_dict, metric)

    def determine_used_prediction_algorithms(self):
        headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
        potential_algorithms = PredictionClass.prediction_methods()
        prediction_algorithms = []
        for algorithm in potential_algorithms:
            if "{} MT Score".format(algorithm) in headers or "{} Score".format(algorithm) in headers:
                prediction_algorithms.append(algorithm)
        return prediction_algorithms

    def determine_columns_used_for_aggregation(self, prediction_algorithms):
        used_columns = [
            "Chromosome", "Start", "Stop", "Reference", "Variant",
            "Transcript", "Variant Type", "Mutation",
            "Protein Position", "Gene Name", "HLA Allele",
            "Mutation Position", "MT Epitope Seq", "WT Epitope Seq",
            "Tumor DNA VAF", "Tumor RNA Depth",
            "Tumor RNA VAF", "Gene Expression", "Transcript Expression",
            "Median MT Score", "Median WT Score", "Median MT Percentile", "Median WT Percentile",
        ]
        for algorithm in prediction_algorithms:
            used_columns.extend(["{} WT Score".format(algorithm), "{} MT Score".format(algorithm), "{} WT Percentile".format(algorithm), "{} MT Percentile".format(algorithm)])
        return used_columns

    def set_column_types(self, prediction_algorithms):
        dtypes = {
            'Chromosome': str,
            "Start": "int32",
            "Stop": "int32",
            'Reference': str,
            'Variant': str,
            "Variant Type": "category",
            "Mutation Position": "category",
            "Median MT Score": "float32",
            "Median MT Percentile": "float16",
            "Protein Position": "str",
        }
        for algorithm in prediction_algorithms:
            if algorithm == 'SMM' or algorithm == 'SMMPMBEC':
                continue
            dtypes["{} MT Score".format(algorithm)] = "float32"
            dtypes["{} MT Percentile".format(algorithm)] = "float16"
        return dtypes

    def execute(self):
        prediction_algorithms = self.determine_used_prediction_algorithms()
        used_columns = self.determine_columns_used_for_aggregation(prediction_algorithms)
        dtypes = self.set_column_types(prediction_algorithms)

        ## get a list of all represented hla types
        hla_types = pd.read_csv(self.input_file, delimiter="\t", usecols=["HLA Allele"])['HLA Allele'].unique()

        ## get a list of unique mutations
        keys = self.get_list_unique_mutation_keys()

        ##do a crude estimate of clonal vaf/purity
        vaf_clonal = self.calculate_clonal_vaf()

        if vaf_clonal is not None:
            metrics = {
                'tumor_purity': self.tumor_purity,
                'vaf_clonal': round(vaf_clonal, 3),
                'vaf_subclonal': round(vaf_clonal/2, 3)
            }
        else:
            metrics = {}

        data = []
        all_epitopes_df = self.read_input_file(used_columns, dtypes)
        for key in keys:
            (df, key_str) = self.get_sub_df(all_epitopes_df, key)
            (best_mut_line, metrics_for_key) = self.get_best_mut_line(df, key_str, hla_types, prediction_algorithms, vaf_clonal, 1000)
            data.append(best_mut_line)
            metrics[key_str] = metrics_for_key
        peptide_table = pd.DataFrame(data=data)
        peptide_table = self.sort_table(peptide_table)

        peptide_table.to_csv(self.output_file, sep='\t', na_rep='NA', index=False, float_format='%.3f')

        self.write_metrics_file(metrics)
        self.copy_pvacview_r_files()



class PvacseqAggregateAllEpitopes(AggregateAllEpitopes, metaclass=ABCMeta):
    def __init__(self, input_file, output_file, tumor_purity=None):
        self.input_file = input_file
        self.output_file = output_file
        self.tumor_purity = tumor_purity
        self.metrics_file = output_file.replace('.tsv', '.metrics.json')

    def get_list_unique_mutation_keys(self):
        key_columns = {
            'Chromosome': str,
            'Start': "int32",
            'Stop': "int32",
            'Reference': str,
            'Variant': str
        }
        key_df = pd.read_csv(self.input_file, delimiter="\t", usecols=key_columns.keys(), dtype=key_columns)
        keys = key_df[['Chromosome', 'Start', 'Stop', 'Reference', 'Variant']].values.tolist()
        keys = [list(i) for i in set(tuple(i) for i in keys)]
        return sorted(keys)

    def calculate_clonal_vaf(self):
        if self.tumor_purity:
            vaf_clonal =  self.tumor_purity * 0.5
            print("Tumor clonal VAF estimated as {} (calculated from user-provided tumor purity of {}). Assuming variants with VAF < {} are subclonal".format(round(vaf_clonal, 3), round(self.tumor_purity, 3), round(vaf_clonal/2, 3)))
            return vaf_clonal
        else:
        #if no tumor purity is provided, make a rough estimate by taking the list of VAFs < 0.6 (assumption is that these are CN-neutral) and return the largest as the marker of the founding clone
            vafs = np.sort(pd.read_csv(self.input_file, delimiter="\t", usecols=["Tumor DNA VAF"])['Tumor DNA VAF'].unique())[::-1]
            vafs_clonal = list(filter(lambda vaf: vaf < 0.6, vafs))
            if len(vafs_clonal) == 0:
                vaf_clonal = 0.6
            else:
                vaf_clonal = vafs_clonal[0]
            print("Tumor clonal VAF estimated as {} (estimated from Tumor DNA VAF data). Assuming variants with VAF < {} are subclonal".format(round(vaf_clonal, 3), round(vaf_clonal/2, 3)))
            return vaf_clonal

    def read_input_file(self, used_columns, dtypes):
        return pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False, na_values="NA", keep_default_na=False, usecols=used_columns, dtype=dtypes)

    def get_sub_df(self, all_epitopes_df, key):
        key_str = "{}-{}-{}-{}-{}".format(key[0], key[1], key[2], key[3], key[4])
        df = (all_epitopes_df[lambda x: (x['Chromosome'] == key[0]) & (x['Start'] == key[1]) & (x['Stop'] == key[2]) & (x['Reference'] == key[3]) & (x['Variant'] == key[4])]).copy()
        df.fillna(value={"Tumor RNA Depth": 0, "Tumor RNA VAF": 0, "Tumor DNA VAF": 0, "Gene Expression": 0}, inplace=True)
        df['Variant Type'] = df['Variant Type'].cat.add_categories('NA')
        df['Mutation Position'] = df['Mutation Position'].cat.add_categories('NA')
        df.fillna(value="NA", inplace=True)
        df['annotation'] = df[['Transcript', 'Gene Name', 'Mutation', 'Protein Position']].agg('-'.join, axis=1)
        df['key'] = key_str
        return (df, key_str)

    def get_best_binder(self, df):
        df.sort_values(by=["Median MT Score", "Median WT Score"], inplace=True, ascending=[True, False])
        return df.iloc[0].to_dict()

    #assign mutations to a "Classification" based on their favorability
    def get_tier(self, mutation, vaf_clonal):
        anchor_residue_pass = True
        anchors = [1, 2, len(mutation["MT Epitope Seq"])-1, len(mutation["MT Epitope Seq"])]
        position = mutation["Mutation Position"]
        if position != "NA":
            if int(float(position)) in anchors:
                if mutation["Median WT Score"] == "NA":
                      anchor_residue_pass = False
                elif mutation["Median WT Score"] < 1000:
                      anchor_residue_pass = False

        #writing these out as explicitly as possible for ease of understanding
        if (mutation["Median MT Score"] < 500 and
           mutation["Tumor RNA VAF"] * mutation["Gene Expression"] > 3 and
           mutation["Tumor DNA VAF"] >= (vaf_clonal/2) and
           anchor_residue_pass):
            return "Pass"

        #relax mt and expr
        if (mutation["Median MT Score"] < 1000 and
           mutation["Tumor RNA VAF"] * mutation["Gene Expression"] > 1 and
           mutation["Tumor DNA VAF"] >= (vaf_clonal/2) and
           anchor_residue_pass):
            return "Relaxed"

        #anchor residues
        if (mutation["Median MT Score"] < 1000 and
           mutation["Tumor RNA VAF"] * mutation["Gene Expression"] > 1 and
           mutation["Tumor DNA VAF"] >= (vaf_clonal/2) and
           not anchor_residue_pass):
            return "Anchor"

        #not in founding clone
        if (mutation["Median MT Score"] < 1000 and
           mutation["Tumor RNA VAF"] * mutation["Gene Expression"] > 1 and
           mutation["Tumor DNA VAF"] < (vaf_clonal/2) and
           anchor_residue_pass):
            return "Subclonal"

        #relax expression.  Include sites that have reasonable vaf but zero overall gene expression
        lowexpr=False
        if ((mutation["Tumor RNA VAF"] * mutation["Gene Expression"] > 0) or
           (mutation["Gene Expression"] == 0 and
           mutation["Tumor RNA Depth"] > 50 and
           mutation["Tumor RNA VAF"] > 0.10)):
             lowexpr=True

        #if low expression is the only strike against it, it gets lowexpr label (multiple strikes will pass through to poor)
        if (mutation["Median MT Score"] < 1000 and
            lowexpr==True and
            mutation["Tumor DNA VAF"] >= (vaf_clonal/2) and
            anchor_residue_pass):
            return "LowExpr"

        #zero expression
        if (mutation["Gene Expression"] == 0 or mutation["Tumor RNA VAF"] == 0) and lowexpr==False:
            return "NoExpr"

        #everything else
        return "Poor"

    def get_good_binders(self, df, max_ic50):
        return df[df["Median MT Score"] < max_ic50]

    def get_unique_good_binders(self, good_binders):
        return pd.DataFrame(good_binders.groupby(['HLA Allele', 'MT Epitope Seq']).size().reset_index())

    def get_good_binders_metrics(self, good_binders, prediction_algorithms, hla_types):
        peptides = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        good_peptides = good_binders["MT Epitope Seq"].unique()
        good_transcripts = good_binders['annotation'].unique()
        for annotation in good_transcripts:
            good_binders_annotation = good_binders[good_binders['annotation'] == annotation]
            for peptide in good_peptides:
                good_binders_peptide_annotation = good_binders_annotation[good_binders_annotation['MT Epitope Seq'] == peptide]
                if len(good_binders_peptide_annotation) > 0:
                    individual_ic50_calls = { 'algorithms': prediction_algorithms }
                    individual_percentile_calls = { 'algorithms': prediction_algorithms }
                    for peptide_type in ['MT', 'WT']:
                        ic50s = {}
                        percentiles = {}
                        ic50_calls = {}
                        percentile_calls = {}
                        for index, line in good_binders_peptide_annotation.to_dict(orient='index').items():
                            ic50s[line['HLA Allele']] = line['Median {} Score'.format(peptide_type)]
                            percentiles[line['HLA Allele']] = line['Median {} Percentile'.format(peptide_type)]
                            ic50_calls[line['HLA Allele']] = [line["{} {} Score".format(algorithm, peptide_type)] for algorithm in prediction_algorithms]
                            percentile_calls[line['HLA Allele']] = [line["{} {} Percentile".format(algorithm, peptide_type)] for algorithm in prediction_algorithms]
                        sorted_ic50s = []
                        sorted_percentiles = []
                        for hla_type in sorted(hla_types):
                            if hla_type in ic50s:
                                sorted_ic50s.append(ic50s[hla_type])
                            else:
                                sorted_ic50s.append('X')
                            if hla_type in percentiles:
                                sorted_percentiles.append(percentiles[hla_type])
                            else:
                                sorted_percentiles.append('X')
                        peptides[annotation][peptide]['ic50s_{}'.format(peptide_type)] = sorted_ic50s
                        peptides[annotation][peptide]['percentiles_{}'.format(peptide_type)] = sorted_percentiles
                        individual_ic50_calls[peptide_type] = ic50_calls
                        individual_percentile_calls[peptide_type] = percentile_calls
                    peptides[annotation][peptide]['hla_types'] = sorted(hla_types)
                    peptides[annotation][peptide]['mutation_position'] = str(good_binders_peptide_annotation.iloc[0]['Mutation Position'])
                    peptides[annotation][peptide]['individual_ic50_calls'] = individual_ic50_calls
                    peptides[annotation][peptide]['individual_percentile_calls'] = individual_percentile_calls
                    wt_peptide = good_binders_peptide_annotation.iloc[0]['WT Epitope Seq']
                    if wt_peptide == 'NA':
                        variant_type = good_binders_peptide_annotation.iloc[0]['Variant Type']
                        if variant_type == 'FS':
                            wt_peptide = 'FS-NA'
                        elif variant_type == 'inframe_ins':
                            wt_peptide = 'INS-NA'
                        elif variant_type == 'inframe_deletion':
                            wt_peptide = 'DEL-NA'
                    peptides[annotation][peptide]['wt_peptide'] = wt_peptide
        return peptides

    def calculate_annotation_count(self, peptides):
        return len(peptides.keys())

    def calculate_unique_peptide_count(self, good_binders):
        return len(good_binders["MT Epitope Seq"].unique())

    def get_default_annotation_count(self):
        return 0

    def get_best_aa_change(self, best):
        if best['Variant Type'] == 'FS':
            return 'FS{}'.format(best['Protein Position'])
        else:
            (wt_aa, mt_aa) = best["Mutation"].split("/")
            return "".join([wt_aa, best["Protein Position"], mt_aa])

    def calculate_allele_expr(self, line):
        if line['Gene Expression'] == 'NA' or line['Tumor RNA VAF'] == 'NA':
            return 'NA'
        else:
            return round(float(line['Gene Expression']) * float(line['Tumor RNA VAF']), 3)

    def assemble_result_line(self, best, key, vaf_clonal, hla, anno_count, peptide_count):
        allele_expr = self.calculate_allele_expr(best)
        tier = self.get_tier(mutation=best, vaf_clonal=vaf_clonal)

        out_dict = { 'ID': key }
        out_dict.update({ k.replace('HLA-', ''):v for k,v in sorted(hla.items()) })
        out_dict.update({
            'Gene': best["Gene Name"],
            'AA Change': self.get_best_aa_change(best),
            'Num Passing Transcripts': anno_count,
            'Best Peptide': best["MT Epitope Seq"],
            'Pos': best["Mutation Position"],
            'Num Passing Peptides': peptide_count,
            'IC50 MT': best["Median MT Score"],
            'IC50 WT': best["Median WT Score"],
            '%ile MT': best["Median MT Percentile"],
            '%ile WT': best["Median WT Percentile"],
            'RNA Expr': best["Gene Expression"],
            'RNA VAF': best["Tumor RNA VAF"],
            'Allele Expr': allele_expr,
            'RNA Depth': best["Tumor RNA Depth"],
            'DNA VAF': best["Tumor DNA VAF"],
            'Tier': tier,
            'Evaluation': 'Pending',
        })
        return out_dict

    def get_metrics(self, df, peptides, best):
        transcripts = list(peptides.keys())
        return {
            'good_binders': peptides,
            'good_binders_transcripts': transcripts,
            'transcript_expr': [df[df["annotation"] == x]['Transcript Expression'].iloc[0] for x in transcripts],
            'DNA VAF': float(best['Tumor DNA VAF']),
            'RNA VAF': float(best['Tumor RNA VAF']),
            'gene_expr': float(best['Gene Expression']),
            'best_peptide_mt': best['MT Epitope Seq'],
            'best_peptide_wt': best['WT Epitope Seq'],
            'best_hla_allele': best['HLA Allele'],
        }

    def write_metrics_file(self, metrics):
        with open(self.metrics_file, 'w') as fh:
            json.dump(metrics, fh, indent=2, separators=(',', ': '))

    #sort the table in our preferred manner
    def sort_table(self, df):
        #make sure the tiers sort in the expected order
        tier_sorter = ["Pass", "Relaxed", "LowExpr", "Anchor", "Subclonal", "Poor", "NoExpr"]
        sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
        df["rank_tier"] = df['Tier'].map(sorter_index)

        df["rank_ic50"] = df["IC50 MT"].rank(ascending=True, method='dense')
        df["rank_expr"] = df["Allele Expr"].rank(ascending=False, method='dense')
        df["rank"] = df["rank_ic50"] + df["rank_expr"]

        df.sort_values(by=["rank_tier", "rank", "Gene", "AA Change"], inplace=True, ascending=True)

        df.drop(labels='rank_tier', axis=1, inplace=True)
        df.drop(labels='rank_ic50', axis=1, inplace=True)
        df.drop(labels='rank_expr', axis=1, inplace=True)
        df.drop(labels='rank', axis=1, inplace=True)

        return df

    def copy_pvacview_r_files(self):
        module_dir = os.path.dirname(__file__)
        r_folder = os.path.abspath(os.path.join(module_dir,"..","tools","pvacview"))
        destination = os.path.abspath(os.path.dirname(self.output_file))
        os.makedirs(os.path.join(destination, "www"), exist_ok=True)
        for i in ["ui.R", "app.R", "server.R", "styling.R", "anchor_and_helper_functions.R"]:
            shutil.copy(os.path.join(r_folder, i), os.path.join(destination, i))
        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            shutil.copy(os.path.join(r_folder, "www", i), os.path.join(destination, "www", i))


class UnmatchedSequenceAggregateAllEpitopes(AggregateAllEpitopes, metaclass=ABCMeta):
    def get_list_unique_mutation_keys(self):
        key_df = pd.read_csv(self.input_file, delimiter="\t", usecols=["Mutation"], dtype={"Mutation": str})
        keys = key_df["Mutation"].values.tolist()
        return sorted(list(set(keys)))

    def calculate_clonal_vaf(self):
        return None

    def read_input_file(self, used_columns, dtypes):
        return pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False, na_values="NA", keep_default_na=False, dtype={"Mutation": str})

    def get_sub_df(self, all_epitopes_df, key):
        df = (all_epitopes_df[lambda x: (x['Mutation'] == key)]).copy()
        return (df, key)

    def get_best_binder(self, df):
        df.sort_values(by=["Median Score"], inplace=True, ascending=True)
        return df.iloc[0]

    def get_tier(self, mutation, vaf_clonal):
        return "NA"

    def get_good_binders(self, df, max_ic50):
        return df[df["Median Score"] < max_ic50]

    def get_unique_good_binders(self, good_binders):
        return pd.DataFrame(good_binders.groupby(['HLA Allele', 'Epitope Seq']).size().reset_index())

    def get_good_binders_metrics(self, good_binders, prediction_algorithms, hla_types):
        return None

    def calculate_annotation_count(self, peptides):
        return "NA"

    def calculate_unique_peptide_count(self, good_binders):
        return len(good_binders["Epitope Seq"].unique())

    def get_default_annotation_count(self):
        return "NA"

    def get_best_aa_change(self, best):
        return 'NA'

    def calculate_allele_expr(self, line):
        return 'NA'

    def assemble_result_line(self, best, key, vaf_clonal, hla, anno_count, peptide_count):
        allele_expr = self.calculate_allele_expr(best)
        tier = self.get_tier(mutation=best, vaf_clonal=vaf_clonal)

        out_dict = { 'ID': key }
        out_dict.update({ k.replace('HLA-', ''):v for k,v in sorted(hla.items()) })
        if 'Gene Name' in best:
            gene = best['Gene Name']
        else:
            gene = 'NA'
        out_dict.update({
            'Gene': gene,
            'AA Change': self.get_best_aa_change(best),
            'Num Passing Transcripts': anno_count,
            'Best Peptide': best["Epitope Seq"],
            'Pos': "NA",
            'Num Passing Peptides': peptide_count,
            'IC50 MT': best["Median Score"],
            'IC50 WT': "NA",
            '%ile MT': best["Median Percentile"],
            '%ile WT': "NA",
            'RNA Expr': "NA",
            'RNA VAF': "NA",
            'Allele Expr': allele_expr,
            'RNA Depth': "NA",
            'DNA VAF': "NA",
            'Tier': tier,
            'Evaluation': 'Pending',
            'ID':key,
        })
        return out_dict

    def get_metrics(self, df, peptides, best):
        return None

    def write_metrics_file(self, metrics):
        pass

    #sort the table in our preferred manner
    def sort_table(self, df):
        df.sort_values(by=["IC50 MT", "ID"], inplace=True, ascending=[True, True])
        return df

    def copy_pvacview_r_files(self):
        pass
