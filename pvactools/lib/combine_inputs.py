import pandas as pd

class CombineInputs():
    def __init__(self, **kwargs):
        self.junctions_df = kwargs['junctions_df']
        self.variants     = kwargs['variant_file']
        self.output_file  = kwargs['output_file']

    def add_junction_coordinates_to_variants(self):
        # read in df
        var_df = pd.read_csv(self.variants, sep='\t')

        # remove version number in annotated to compare with filtered junctions file
        var_df[['transcript_id', 'transcript_version']] = var_df['transcript_name'].str.split('.', expand=True)
        var_df = var_df.astype({'transcript_version': 'float64'})

        # # create new cols
        # var_df[['junction_variant_start', 'junction_variant_stop']] = 0

        # # set up variant_category
        # var_df['variant_category'] = 'SNV'
        # var_df.loc[var_df['reference'].str.len() > var_df['variant'].str.len(), 'variant_category'] = 'DEL'
        # var_df.loc[var_df['reference'].str.len() < var_df['variant'].str.len(), 'variant_category'] = 'INS'

        # # copy values - SNVs
        # var_df.loc[var_df['variant_category'] == 'SNV', 'junction_variant_start'] = var_df['start']
        # var_df.loc[var_df['variant_category'] == 'SNV', 'junction_variant_stop'] = var_df['stop']

        # # deletions - bp size matters
        # var_df.loc[var_df['variant_category'] == 'DEL', 'junction_variant_start'] = var_df['start']
        # var_df.loc[var_df['variant_category'] == 'DEL', 'junction_variant_stop'] = var_df['start'] + var_df['reference'].str.len() - var_df['variant'].str.len()

        # # insertions
        # var_df.loc[var_df['variant_category'] == 'INS', 'junction_variant_start'] = var_df['start']
        # var_df.loc[var_df['variant_category'] == 'INS', 'junction_variant_stop'] = var_df['start']
        
        # MNV support (exploded variants now in SNV notation)
        
        var_df = var_df.rename(columns={'ensembl_gene_id': 'gene_id'}).drop(columns=['transcript_support_level'])

        # format junction variant info to match vcf
        var_df['variant_info'] = var_df['chromosome_name'] + ':' + var_df['start'].astype('string') + '-' + var_df['stop'].astype('string')
        
        return var_df

    def merge_and_write(self, j_df, var_df):
        # merge cols: 'transcript_support_level', 'gene_name', 'start' ('stop' vs. 'end')
        # merge by transcript and variant coors
        merged_df = j_df.merge(var_df, on=['transcript_id', 'transcript_version', 'gene_name', 'gene_id', 'variant_info']).drop_duplicates()

        # create index to match with kmers
        merged_df['index'] = merged_df['gene_name'] + '.' + merged_df['transcript_id'] + '.' + merged_df['name'] + '.' + merged_df['variant'] + '.' + merged_df['anchor']
        
        # cols for frameshift info
        merged_df[['wt_protein_length', 'alt_protein_length', 'frameshift_event']] = 'NA'
        
        # create final file
        merged_df.to_csv(self.output_file, sep='\t', index=False)

        return merged_df

    def execute(self):
        # create dfs
        variant_df = self.add_junction_coordinates_to_variants()
        # merge dfs and create associated combined file
        combined_df = self.merge_and_write(self.junctions_df, variant_df)

        # is protein is NA, can I go ahead and remove the lines? do these correspond to skipped junctions in jtf?
        
        return combined_df

