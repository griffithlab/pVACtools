import os
import pandas as pd

class CombineInputs():
    def __init__(self, **kwargs):
        self.junctions   = kwargs['junctions_file']
        self.variants    = kwargs['variant_file']
        self.sample_name = kwargs['sample_name']
        self.output_dir  = kwargs['output_dir']
        self.output_prefix = f'{self.output_dir}/{self.sample_name}'

    def add_junction_coordinates_to_variants(self):
        # read in df
        var_df = pd.read_csv(self.variants, sep='\t')

        # remove version number in annotated to compare with filtered junctions file
        var_df['transcript_name'] = var_df['transcript_name'].str.split('.', expand=True)[[0]]

        # create new cols
        var_df[['junction_variant_start', 'junction_variant_stop']] = 0

        # set up variant_category
        var_df['variant_category'] = 'SNV'
        var_df.loc[var_df['reference'].str.len() > var_df['variant'].str.len(), 'variant_category'] = 'DEL'
        var_df.loc[var_df['reference'].str.len() < var_df['variant'].str.len(), 'variant_category'] = 'INS'

        # copy values - SNVs
        var_df.loc[var_df['variant_category'] == 'SNV', 'junction_variant_start'] = var_df['start']
        var_df.loc[var_df['variant_category'] == 'SNV', 'junction_variant_stop'] = var_df['stop']

        # deletions
        var_df.loc[var_df['variant_category'] == 'DEL', 'junction_variant_start'] = var_df['start']
        var_df.loc[var_df['variant_category'] == 'DEL', 'junction_variant_stop'] = var_df['start'] + 1

        # insertions
        var_df.loc[var_df['variant_category'] == 'INS', 'junction_variant_start'] = var_df['start'] -1
        var_df.loc[var_df['variant_category'] == 'INS', 'junction_variant_stop'] = var_df['start']

        # MNVs
        var_df.loc[(var_df['variant_category'] == 'MNV'), 'junction_variant_start'] = var_df['start']
        var_df.loc[var_df['variant_category'] == 'MNV', 'junction_variant_stop'] = var_df['stop']

        # format to match junctions
        var_df['variant_info'] = var_df['chromosome_name'] + ':' + var_df['junction_variant_start'].astype('string') + '-' + var_df['junction_variant_stop'].astype('string')

        var_df.to_csv(f'{self.output_prefix}_annotated.tsv', sep='\t', index=False)
        
        return var_df

    def execute(self):

        variant_df = self.add_junction_coordinates_to_variants()
        junction_df = pd.read_csv(f'{self.output_prefix}_filtered.tsv', sep='\t')

        # also on transcripts
        merged_df = junction_df.merge(variant_df, on=['variant_info']).drop_duplicates()

        merged_df.to_csv(f'{self.output_prefix}_combined.tsv', sep='\t', index=False)

# variant_df unique variants = 46
# junction_df unique variants = 13
# merged_df = 12 - missing MNV