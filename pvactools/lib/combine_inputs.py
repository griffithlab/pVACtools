import sys
import pandas as pd


class CombineInputs:
    def __init__(self, **kwargs):
        self.junctions_df = kwargs['junctions_df']
        self.variants     = kwargs['variant_file']
        self.output_file  = kwargs['output_file']
        self.output_dir   = kwargs['output_dir']

    def add_junction_coordinates_to_variants(self):
        # read in df
        var_df = pd.read_csv(self.variants, sep='\t')

        # remove version number in annotated to compare with filtered junctions file
        var_df[['transcript_id', 'transcript_version']] = var_df['transcript_name'].str.split('.', expand=True)
        var_df = var_df.loc[var_df['transcript_id'].str.startswith('ENST') == True]
        var_df['transcript_version'] = var_df['transcript_version'].astype('int64')

        # create new cols
        var_df[['variant_start', 'variant_stop']] = 0

        # set up variant_category
        var_df['variant_category'] = 'SNV'
        var_df.loc[var_df['reference'].str.len() > var_df['variant'].str.len(), 'variant_category'] = 'DEL'
        var_df.loc[var_df['reference'].str.len() < var_df['variant'].str.len(), 'variant_category'] = 'INS'

        # copy values - SNVs
        var_df.loc[var_df['variant_category'] == 'SNV', 'variant_start'] = var_df['start']
        var_df.loc[var_df['variant_category'] == 'SNV', 'variant_stop'] = var_df['stop']

        # deletions - bp size matters
        var_df.loc[var_df['variant_category'] == 'DEL', 'variant_start'] = var_df['start'] - 1
        var_df.loc[var_df['variant_category'] == 'DEL', 'variant_stop'] = var_df['start']

        # insertions
        var_df.loc[var_df['variant_category'] == 'INS', 'variant_start'] = var_df['start']
        var_df.loc[var_df['variant_category'] == 'INS', 'variant_stop'] = var_df['start']

        # MNV support (exploded variants now in SNV notation)

        var_df = var_df.rename(columns={'ensembl_gene_id': 'gene_id'}).drop(columns=['transcript_support_level'])

        # format junction variant info to match vcf
        var_df['variant_info'] = var_df['chromosome_name'] + ':' + var_df['variant_start'].astype('string') + '-' + var_df['variant_stop'].astype('string')

        return var_df

    def merge_and_write(self, j_df, var_df):
        # if protein change/seq is NA in var_df, go ahead and remove the lines bc if there is no protein change, then can't create alt transcript
        j_df['transcript_version'] = j_df['transcript_version'].astype('int64')

        # removed gene_name, gene_id - do these need to be skipped?
        merged_df = j_df.merge(var_df, on=['gene_name', 'transcript_id', 'transcript_version', 'variant_info', 'gene_id'])

        # check that everything is merging
        left_merge = j_df.merge(var_df, on=['variant_info', 'gene_name', 'transcript_id', 'transcript_version', 'gene_id'], how='left', indicator=True)
        not_merged_lines = left_merge.loc[left_merge['_merge'] != 'both']
        if not not_merged_lines.empty:
            # warning: if there are any that don't merge,
            print(not_merged_lines[['variant_info', 'gene_name', 'transcript_id', 'transcript_version']])
            # make this back into a warning
            print(
                'Warning: The above variant/transcript/gene combination is present in the junctions file, but not in the VEP-annotated VCF. Skipping.'
            )

        # create index to match with kmers
        merged_df['index'] = merged_df['gene_name'] + '.' + merged_df['transcript_id'] + '.' + merged_df['name'] + '.' + merged_df['variant_info'] + '.' + merged_df['anchor']

        # cols for frameshift info
        merged_df[['wt_protein_length', 'alt_protein_length', 'frameshift_event']] = pd.NA

        return merged_df

    def execute(self):
        # create dfs
        variant_df = self.add_junction_coordinates_to_variants()
        # merge dfs and create associated combined file
        combined_df = self.merge_and_write(self.junctions_df, variant_df)

        return combined_df
