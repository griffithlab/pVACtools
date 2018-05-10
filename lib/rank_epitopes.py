import pandas as pd

class RankEpitopes:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def headers(self):
        return [
            'Gene Name',
            'Mutation',
            'Protein Position',
            'HGVSc',
            'HGVSp',
            'HLA Allele',
            'MT Epitope Seq',
            'MT IC50',
            'WT IC50',
            'Fold Change',
            'Tumor DNA Depth',
            'Tumor DNA VAF',
            'Tumor RNA Depth',
            'Tumor RNA VAF',
            'Gene Expression',
            'Score',
        ]

    def execute(self):
        df = pd.read_csv(self.input_file, sep='\t', index_col=False)
        df['mt_score_rank'] = df['MT IC50'].rank(numeric_only=True, ascending=False, method='dense').fillna(value=0.0)
        df['fold_change_rank'] = df['Fold Change'].rank(numeric_only=True, ascending=True, method='dense').fillna(value=0.0)
        df['mt_allele_exp'] = df['Tumor RNA VAF'] * df['Gene Expression']
        df['mt_allele_exp_rank'] = df['mt_allele_exp'].rank(numeric_only=True, ascending=True, method='dense').fillna(value=0.0)
        df['tumor_dna_vaf_rank'] = df['Tumor DNA VAF'].rank(numeric_only=True, ascending=True, method='dense').fillna(value=0.0)
        df['Score'] = df['mt_score_rank'] + df['fold_change_rank'] + (df['mt_allele_exp_rank'] * 2) + (df['tumor_dna_vaf_rank'] / 2)
        df.sort_values(by=['Score'], inplace=True, ascending=False)
        df.to_csv(self.output_file, sep='\t', na_rep='NA', columns=self.headers(), index=False)
