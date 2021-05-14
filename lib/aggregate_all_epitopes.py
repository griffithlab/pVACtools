import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import json
from Bio import SeqIO
from .prediction_class import *

class AggregateAllEpitopes:
    def __init__(self, input_file, output_file, file_type='pVACseq', fasta_file=None):
        self.input_file = input_file
        self.output_file = output_file
        self.metrics_file = output_file.replace('.tsv', '.metrics.json')
        self.fasta_file = fasta_file
        self.file_type = file_type

    #assign mutations to a "Classification" based on their favorability
    def get_tier(self, mutation, vaf_clonal):
        anchor_residue_pass = True
        anchors = [1, 2, len(mutation["MT Epitope Seq"])-1, len(mutation["MT Epitope Seq"])]
        if mutation["Mutation Position"] in anchors:
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

    def get_best_mut_line(self, df, hla_types, prediction_algorithms, vaf_clonal, max_ic50=1000):
        #order by best median score and get best ic50 peptide
        if self.file_type == 'pVACbind':
            df.sort_values(by=["Median Score"], inplace=True, ascending=True)
        else:
            df.sort_values(by=["Median MT Score", "Median WT Score"], inplace=True, ascending=[True, False])
        best = df.iloc[0]

        if self.file_type == 'pVACbind':
            tier = "NA"
        else:
            tier = self.get_tier(best, vaf_clonal)

        #these counts should represent only the "good binders" with ic50 < max
        #for all sites other than tier4 slop
        if self.file_type == 'pVACbind':
            good_binders = df[df["Median Score"] < max_ic50]
        else:
            good_binders = df[df["Median MT Score"] < max_ic50]
        if len(good_binders) > 0:
            good_binders_uniq = pd.DataFrame(good_binders.groupby(['HLA Allele', 'MT Epitope Seq']).size().reset_index())
            good_binders_hla = Counter(good_binders_uniq["HLA Allele"])
            hla = dict(map(lambda x : (x, good_binders_hla[x]) if x in good_binders_hla else (x, ""), hla_types))
            #get a list of all unique gene/transcript/aa_change combinations
            #store a count of all unique peptides that passed
            if self.file_type == 'pVACbind':
                anno_count = "NA"
                peptides = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict)))
                for index, line in good_binders.to_dict(orient='index').items():
                    peptides[line['annotation']][line['Epitope Seq']]['ic_50']["{} MT".format(line['HLA Allele'])] = line['Median Score']
                    peptides[line['annotation']][line['Epitope Seq']]['ic_50']["{} WT".format(line['HLA Allele'])] = "NA"
                peptide_count = len(good_binders["Epitope Seq"].unique())
            else:
                peptides = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
                for peptide in good_binders["MT Epitope Seq"].unique():
                    good_binders_peptide = good_binders[good_binders['MT Epitope Seq'] == peptide]
                    individual_ic50_calls = { 'algorithms': prediction_algorithms }
                    individual_percentile_calls = { 'algorithms': prediction_algorithms }
                    for peptide_type in ['MT', 'WT']:
                        ic50s = {}
                        percentiles = {}
                        ic50_calls = {}
                        percentile_calls = {}
                        for index, line in good_binders_peptide.to_dict(orient='index').items():
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
                        peptides[line['annotation']][peptide]['ic50s_{}'.format(peptide_type)] = sorted_ic50s
                        peptides[line['annotation']][peptide]['percentiles_{}'.format(peptide_type)] = sorted_percentiles
                        individual_ic50_calls[peptide_type] = ic50_calls
                        individual_percentile_calls[peptide_type] = percentile_calls
                    peptides[line['annotation']][peptide]['hla_types'] = sorted(hla_types)
                    peptides[line['annotation']][peptide]['mutation_position'] = str(good_binders_peptide.iloc[0]['Mutation Position'])
                    peptides[line['annotation']][peptide]['individual_ic50_calls'] = individual_ic50_calls
                    peptides[line['annotation']][peptide]['individual_percentile_calls'] = individual_percentile_calls
                    wt_peptide = good_binders_peptide.iloc[0]['WT Epitope Seq']
                    if wt_peptide == 'NA':
                        variant_type = good_binders_peptide.iloc[0]['Variant Type']
                        if variant_type == 'FS':
                            wt_peptide = 'FS-NA'
                        elif variant_type == 'inframe_ins':
                            wt_peptide = 'INS-NA'
                        elif variant_type == 'inframe_deletion':
                            wt_peptide = 'DEL-NA'
                    peptides[line['annotation']][peptide]['wt_peptide'] = wt_peptide
                anno_count = len(peptides.keys())
                peptide_count = len(good_binders["MT Epitope Seq"].unique())
        else:
            hla = dict(map(lambda x : (x, ""), hla_types))
            peptides = {}
            anno_count = 0
            peptide_count = 0

        if self.file_type == 'bedpe':
            best['aachange'] = best['key']
        elif self.file_type == 'pVACbind':
            best['aachange'] = best['Mutation']
        else:
            if best['Variant Type'] == 'FS':
                best['aachange'] = 'FS{}'.format(best['Protein Position'])
            else:
                (wt_aa, mt_aa) = best["Mutation"].split("/")
                best["aachange"] = "".join([wt_aa, best["Protein Position"], mt_aa])

        allele_expr = self.calculate_allele_expr(best)

        #assemble the line
        out_dict = { k.replace('HLA-', ''):v for k,v in hla.items() }
        if self.file_type == 'pVACbind':
            out_dict.update({
                'Gene': ["NA"],
                'AA Change': [best["aachange"]],
                'Num Passing Transcripts': [anno_count],
                'Best Peptide': [best["Epitope Seq"]],
                'Pos': ["NA"],
                'Num Passing Peptides': [peptide_count],
                'IC50 MT': [best["Median Score"]],
                'IC50 WT': ["NA"],
                '%ile MT': [best["Median Percentile"]],
                '%ile WT': ["NA"],
                'RNA Expr': ["NA"],
                'RNA VAF': ["NA"],
                'Allele Expr': allele_expr,
                'RNA Depth': ["NA"],
                'DNA VAF': ["NA"],
                'Tier': [tier],
                'Evaluation': 'Pending',
                'ID':[best['key']],
            })
        else:
            out_dict.update({
                'Gene': [best["Gene Name"]],
                'AA Change': [best["aachange"]],
                'Num Passing Transcripts': [anno_count],
                'Best Peptide': [best["MT Epitope Seq"]],
                'Pos': [best["Mutation Position"]],
                'Num Passing Peptides': [peptide_count],
                'IC50 MT': [best["Median MT Score"]],
                'IC50 WT': [best["Median WT Score"]],
                '%ile MT': [best["Median MT Percentile"]],
                '%ile WT': [best["Median WT Percentile"]],
                'RNA Expr': [best["Gene Expression"]],
                'RNA VAF': [best["Tumor RNA VAF"]],
                'Allele Expr': allele_expr,
                'RNA Depth': [best["Tumor RNA Depth"]],
                'DNA VAF': [best["Tumor DNA VAF"]],
                'Tier': [tier],
                'Evaluation': 'Pending',
                'ID':[best['key']],
            })

        df_out = pd.DataFrame.from_dict(out_dict)
        return (df_out, peptides)

    def calculate_allele_expr(self, line):
        if self.file_type == 'pVACbind':
            return 'NA'
        elif line['Gene Expression'] == 'NA' or line['Tumor RNA VAF'] == 'NA':
            return 'NA'
        else:
            return round(float(line['Gene Expression']) * float(line['Tumor RNA VAF']), 3)

    #sort the table in our preferred manner
    def sort_table(self, df):
        if self.file_type == 'pVACbind':
            df.sort_values(by=["IC50 MT"], inplace=True, ascending=True)
        else:
            #make sure the tiers sort in the expected order
            tier_sorter = ["Pass", "Relaxed", "LowExpr", "Anchor", "Subclonal", "Poor", "NoExpr"]
            sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
            df["rank_tier"] = df['Tier'].map(sorter_index)

            df["rank_ic50"] = df["IC50 MT"].rank(ascending=True, method='dense')
            df["expr"] = df["RNA Expr"] * df["RNA VAF"]
            df["rank_expr"] = df["expr"].rank(ascending=False, method='dense')
            df["rank"] = df["rank_ic50"] + df["rank_expr"]

            df.sort_values(by=["rank_tier", "rank", "Gene", "AA Change"], inplace=True, ascending=True)

            df.drop('rank_tier', 1, inplace=True)
            df.drop('rank_ic50', 1, inplace=True)
            df.drop('expr', 1, inplace=True)
            df.drop('rank_expr', 1, inplace=True)
            df.drop('rank', 1, inplace=True)

        return df

    def parse_fasta_file(self):
        seq_dict = {}
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            seq_dict[record.id] = str(record.seq)
        return seq_dict

    def execute(self):
        peptide_fastas = self.parse_fasta_file()
        df = pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False)
        df.fillna(value={"Tumor RNA Depth": 0, "Tumor RNA VAF": 0, "Tumor DNA VAF": 0, "Gene Expression": 0}, inplace=True)
        df.fillna(value="NA", inplace=True)

        ## get a list of all prediction algorithms
        potential_algorithms = PredictionClass.prediction_methods()
        headers = list(df.columns)
        prediction_algorithms = []
        for algorithm in potential_algorithms:
            if "{} MT Score".format(algorithm) in headers or "{} Score".format(algorithm) in headers:
                prediction_algorithms.append(algorithm)


        ## get a list of all represented hla types
        hla_types = df['HLA Allele'].unique()

        ## get a list of unique mutations
        if self.file_type == 'pVACbind':
            df["key"] = df["Mutation"]
            df['annotation'] = df['Mutation']
            vaf_clonal = None
        else:
            for column in ['Chromosome', 'Start', 'Stop', 'Protein Position', 'Mutation']:
                df[column] = df[column].astype(str)
            if self.file_type == 'bedpe':
                df["key"] = df[['Chromosome', 'Start', 'Stop']].agg(' | '.join, axis=1)
            else:
                df["key"] = df[['Chromosome', 'Start', 'Stop', 'Reference', 'Variant']].agg('-'.join, axis=1)
            df['annotation'] = df[['Transcript', 'Gene Name', 'Mutation', 'Protein Position']].agg('-'.join, axis=1)
            #do a crude estimate of clonal vaf/purity
            vafs = np.sort(df['Tumor DNA VAF'].unique())[::-1]
            vaf_clonal = list(filter(lambda vaf: vaf < 0.6, vafs))[0]

        keys = df["key"].unique()

        columns = ['ID']
        columns.extend(sorted([x.replace('HLA-', '') for x in hla_types]))
        columns.extend(['Gene', 'AA Change', 'Num Passing Transcripts', 'Best Peptide', 'Pos', 'Num Passing Peptides', 'IC50 MT', 'IC50 WT', '%ile MT', '%ile WT', 'RNA Expr', 'RNA VAF', 'Allele Expr', 'RNA Depth', 'DNA VAF', 'Tier', 'Evaluation'])
        peptide_table = pd.DataFrame(columns=columns)
        metrics = {}
        for key in keys:
            df_subset = df[df["key"] == key]
            (best_mut_line, peptides) = self.get_best_mut_line(df_subset, hla_types, prediction_algorithms, vaf_clonal, 1000)
            peptide_table = peptide_table.append(best_mut_line, sort=False)
            all_peptides = defaultdict(lambda: defaultdict(list))
            for index, line in df_subset.to_dict(orient='index').items():
                if self.file_type == 'pVACbind':
                    all_peptides[line['annotation']][line['Epitope Seq']].append(line)
                else:
                    all_peptides[line['annotation']][line['MT Epitope Seq']].append(line)
            transcripts = list(peptides.keys())
            metrics[key] = {
                'good_binders': peptides,
                'good_binders_transcripts': transcripts,
                'transcript_expr': [df_subset[df_subset["annotation"] == x]['Transcript Expression'].iloc[0] for x in transcripts],
                'DNA VAF': float(best_mut_line['DNA VAF']),
                'RNA VAF': float(best_mut_line['RNA VAF']),
                'gene_expr': float(best_mut_line['RNA Expr']),
            }
            if self.fasta_file is not None:
                if self.file_type == 'pVACbind':
                    metrics[key]['peptide'] = peptides[line['Index']]
                else:
                    if line['Variant Type'] == 'FS':
                        index = '%s.%s.%s.%s%s/%s' % (line['Gene Name'], line['Transcript'], line['Variant Type'], line['Protein Position'], line['Reference'], line['Variant'])
                    else:
                        index = '%s.%s.%s.%s%s' % (line['Gene Name'], line['Transcript'], line['Variant Type'], line['Protein Position'], line['Mutation'])
                    wt_matches = [i for i in peptide_fastas.keys() if i.endswith(index) and i.startswith('WT')]
                    mt_matches = [i for i in peptide_fastas.keys() if i.endswith(index) and i.startswith('MT')]
                    if len(wt_matches) == 1 and len(mt_matches) == 1:
                        metrics[key]['mt_peptide'] = peptide_fastas[mt_matches[0]]
                        metrics[key]['wt_peptide'] = peptide_fastas[wt_matches[0]]
        peptide_table = self.sort_table(peptide_table)

        peptide_table.to_csv(self.output_file, sep='\t', na_rep='NA', index=False)

        with open(self.metrics_file, 'w') as fh:
            json.dump(metrics, fh, indent=2, separators=(',', ': '))
