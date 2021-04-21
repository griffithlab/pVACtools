import pandas as pd
import numpy as np

class AggregateAllEpitopes:
    def __init__(self, input_file, output_file, file_type='pVACseq'):
        self.input_file = input_file
        self.output_file = output_file
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

    def get_best_mut_line(self, df, hla_types, vaf_clonal, max_ic50=1000):
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
            good_binders_hla = good_binders["HLA Allele"].unique()
            hla = dict(map(lambda x : (x, 'X') if x in good_binders_hla else (x, ""), hla_types))
            #get a list of all unique gene/transcript/aa_change combinations
            #store a count of all unique peptides that passed
            if self.file_type == 'pVACbind':
                anno_count = "NA"
                peptide_count = len(good_binders["Epitope Seq"].unique())
            else:
                anno_count = len(good_binders[['Transcript', 'Gene Name', 'Mutation', 'Protein Position']].agg('-'.join, axis=1).unique())
                peptide_count = len(good_binders["MT Epitope Seq"].unique())
        else:
            hla = dict(map(lambda x : (x, ""), hla_types))
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

        #assemble the line
        out_dict = hla
        if self.file_type == 'pVACbind':
            out_dict.update({
                'Gene': ["NA"],
                'AA_change': [best["aachange"]],
                'Num_Transcript': [anno_count],
                'Peptide': [best["Epitope Seq"]],
                'Pos': ["NA"],
                'Num_Peptides': [peptide_count],
                'ic50_MT': [best["Median Score"]],
                'ic50_WT': ["NA"],
                'percentile_MT': [best["Median Percentile"]],
                'percentile_WT': ["NA"],
                'RNA_expr': ["NA"],
                'RNA_VAF': ["NA"],
                'RNA_Depth': ["NA"],
                'DNA_VAF': ["NA"],
                'tier': [tier],
            })
        else:
            out_dict.update({
                'Gene': [best["Gene Name"]],
                'AA_change': [best["aachange"]],
                'Num_Transcript': [anno_count],
                'Peptide': [best["MT Epitope Seq"]],
                'Pos': [best["Mutation Position"]],
                'Num_Peptides': [peptide_count],
                'ic50_MT': [best["Median MT Score"]],
                'ic50_WT': [best["Median WT Score"]],
                'percentile_MT': [best["Median MT Percentile"]],
                'percentile_WT': [best["Median WT Percentile"]],
                'RNA_expr': [best["Gene Expression"]],
                'RNA_VAF': [best["Tumor RNA VAF"]],
                'RNA_Depth': [best["Tumor RNA Depth"]],
                'DNA_VAF': [best["Tumor DNA VAF"]],
                'tier': [tier],
            })

        df_out = pd.DataFrame.from_dict(out_dict)
        return df_out

    #sort the table in our preferred manner
    def sort_table(self, df):
        if self.file_type == 'pVACbind':
            df.sort_values(by=["ic50_MT"], inplace=True, ascending=True)
        else:
            #make sure the tiers sort in the expected order
            tier_sorter = ["Pass", "Relaxed", "LowExpr", "Anchor", "Subclonal", "Poor", "NoExpr"]
            sorter_index = dict(zip(tier_sorter,range(len(tier_sorter))))
            df["rank_tier"] = df['tier'].map(sorter_index)

            df["rank_ic50"] = df["ic50_MT"].rank(ascending=True, method='dense')
            df["expr"] = df["RNA_expr"] * df["RNA_VAF"]
            df["rank_expr"] = df["expr"].rank(ascending=False, method='dense')
            df["rank"] = df["rank_ic50"] + df["rank_expr"]

            df.sort_values(by=["rank_tier", "rank", "Gene", "AA_change"], inplace=True, ascending=True)

            df.drop('rank_tier', 1, inplace=True)
            df.drop('rank_ic50', 1, inplace=True)
            df.drop('expr', 1, inplace=True)
            df.drop('rank_expr', 1, inplace=True)
            df.drop('rank', 1, inplace=True)

        return df

    def execute(self):
        df = pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False, na_values="NA", keep_default_na=False)
        df.fillna(value={"Tumor RNA Depth": 0, "Tumor RNA VAF": 0, "Tumor DNA VAF": 0, "Gene Expression": 0}, inplace=True)

        ## get a list of all represented hla types
        hla_types = df['HLA Allele'].unique()

        ## get a list of unique mutations
        if self.file_type == 'pVACbind':
            df["key"] = df["Mutation"]
            vaf_clonal = None
        else:
            for column in ['Chromosome', 'Start', 'Stop', 'Protein Position', 'Mutation']:
                df[column] = df[column].astype(str)
            if self.file_type == 'bedpe':
                df["key"] = df[['Chromosome', 'Start', 'Stop']].agg(' | '.join, axis=1)
            else:
                df["key"] = df[['Chromosome', 'Start', 'Stop', 'Reference', 'Variant']].agg('-'.join, axis=1)
            #do a crude estimate of clonal vaf/purity
            vafs = np.sort(df['Tumor DNA VAF'].unique())[::-1]
            vaf_clonal = list(filter(lambda vaf: vaf < 0.6, vafs))[0]

        keys = df["key"].unique()

        columns = hla_types.tolist()
        columns.extend(['Gene', 'AA_change', 'Num_Transcript', 'Peptide', 'Pos', 'Num_Peptides', 'ic50_MT', 'ic50_WT', 'percentile_MT', 'percentile_WT', 'RNA_expr', 'RNA_VAF', 'RNA_Depth', 'DNA_VAF', 'tier'])
        peptide_table = pd.DataFrame(columns=columns)
        for key in keys:
            df_subset = df[df["key"] == key]
            best_mut_line = self.get_best_mut_line(df_subset, hla_types, vaf_clonal, 1000)
            peptide_table = peptide_table.append(best_mut_line, sort=False)
        peptide_table = self.sort_table(peptide_table)

        peptide_table.to_csv(self.output_file, sep='\t', na_rep='NA', index=False)
