import pandas as pd

class RankEpitopes:
    def __init__(self, input_file, output_file, top_score_metric):
        self.input_file = input_file
        self.output_file = output_file
        self.top_score_metric = top_score_metric

    def headers(self):
        return [
            'Gene Name',
            'Mutation',
            'Protein Position',
            'HGVSc',
            'HGVSp',
            'HLA Allele',
            'Mutation Position',
            'MT Epitope Seq',
            'Median MT Score',
            'Median WT Score',
            'Median Fold Change',
            'Best MT Score',
            'Corresponding WT Score',
            'Corresponding Fold Change',
            'Tumor DNA Depth',
            'Tumor DNA VAF',
            'Tumor RNA Depth',
            'Tumor RNA VAF',
            'Gene Expression',
            'Rank',
        ]

    def execute(self):
        if self.top_score_metric == 'median':
            score_column = 'Median MT Score'
            fold_change_column = 'Median Fold Change'
        elif self.top_score_metric == 'lowest':
            score_column = 'Best MT Score'
            fold_change_column = 'Corresponding Fold Change'
        df = pd.read_csv(self.input_file, sep='\t', index_col=False)
        df['mt_score_rank'] = df[score_column].rank(numeric_only=True, ascending=False, method='dense').fillna(value=0.0)
        df['fold_change_rank'] = df[fold_change_column].rank(numeric_only=True, ascending=True, method='dense').fillna(value=0.0)
        df['mt_allele_exp'] = df['Tumor RNA VAF'] * df['Gene Expression']
        df['mt_allele_exp_rank'] = df['mt_allele_exp'].rank(numeric_only=True, ascending=True, method='dense').fillna(value=0.0)
        df['tumor_dna_vaf_rank'] = df['Tumor DNA VAF'].rank(numeric_only=True, ascending=True, method='dense').fillna(value=0.0)
        df['score'] = df['mt_score_rank'] + df['fold_change_rank'] + (df['mt_allele_exp_rank'] * 2) + (df['tumor_dna_vaf_rank'] / 2)
        df['Rank'] = df['score'].rank(ascending=False, method='dense').astype(int)
        df.sort_values(by=['Rank'], inplace=True, ascending=True)
        df.to_csv(self.output_file, sep='\t', na_rep='NA', columns=self.headers(), index=False)
