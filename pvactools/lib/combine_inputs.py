import pandas as pd

class CombineInputs():
    def __init__(self, **kwargs):
        self.junctions_df = kwargs['junctions_df']
        self.variants     = kwargs['variant_file']
        #self.output_file  = kwargs['output_file']
        self.maximum_transcript_support_level = kwargs['maximum_transcript_support_level']

    def add_junction_coordinates_to_variants(self):
        # read in df
        var_df = pd.read_csv(self.variants, sep='\t')

        # filter variants by tsl and biotype
        var_df = var_df[(var_df['transcript_support_level'] <= self.maximum_transcript_support_level) & (var_df['biotype'] == 'protein_coding')]

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
        var_df.loc[var_df['variant_category'] == 'DEL', 'junction_variant_start'] = var_df['start'] -1
        var_df.loc[var_df['variant_category'] == 'DEL', 'junction_variant_stop'] = var_df['start']

        # insertions
        var_df.loc[var_df['variant_category'] == 'INS', 'junction_variant_start'] = var_df['start'] -1
        var_df.loc[var_df['variant_category'] == 'INS', 'junction_variant_stop'] = var_df['start']

        # MNV support (exploded variants now in SNV notation)

        # format to match junctions
        var_df['variant_info'] = var_df['chromosome_name'] + ':' + var_df['junction_variant_start'].astype('string') + '-' + var_df['junction_variant_stop'].astype('string')
        
        return var_df

    def merge_and_write(self, j_df, var_df):
        # merge by transcript and variant coors
        merged_df = j_df.merge(var_df, on=['transcript_name', 'variant_info']).drop_duplicates()

        # create index to match with kmers
        merged_df['index'] = merged_df['Gene_name'] + '.' + merged_df['transcript_name'] + '.' + merged_df['name'] + '.' + merged_df['variant_info'] + '.' + merged_df['anchor']
        
        # cols for frameshift info
        merged_df[['wt_protein_length', 'alt_protein_length', 'frameshift_event']] = 'NA'
        
        # create final file
        #merged_df.to_csv(self.output_file, sep='\t', index=False)

        return merged_df

    def execute(self):
        # create dfs
        variant_df = self.add_junction_coordinates_to_variants()
        # merge dfs and create associated combined file
        combined_df = self.merge_and_write(self.junctions_df, variant_df)
        
        return combined_df
