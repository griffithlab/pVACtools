import os
import pandas as pd

# class CombineInputs():
#     def __init__(self, **kwargs):
#         self.junctions   = kwargs['junctions_file']
#         self.variants    = kwargs['variant_file']
#         self.sample_name = kwargs['sample_name']
#         self.output_dir  = kwargs['output_dir']
#         self.output_path = f'{self.output_dir}/{self.sample_name}_combined.tsv'

def add_junction_coordinates_to_variants(variants):
    # read in df
    var_df = pd.read_csv(variants, sep='\t') #self.

    # create new cols
    var_df[['junction_variant_start', 'junction_variant_stop']] = 0

    # set up variant_category
    var_df['variant_category'] = 'SNV'
    var_df.loc[var_df['reference'].str.len() > var_df['variant'].str.len(), 'variant_category'] = 'DEL'
    var_df.loc[var_df['reference'].str.len() < var_df['variant'].str.len(), 'variant_category'] = 'INS'
    # add MNV

    # copy values - SNVs
    var_df.loc[var_df['variant_category'] == 'SNV', 'junction_variant_start'] = var_df['start']
    var_df.loc[var_df['variant_category'] == 'SNV', 'junction_variant_stop'] = var_df['stop']

    # deletions
    var_df.loc[var_df['variant_category'] == 'DEL', 'junction_variant_start'] = var_df['start']
    var_df.loc[var_df['variant_category'] == 'DEL', 'junction_variant_stop'] = var_df['start'] + 1

    # insertions
    var_df.loc[var_df['variant_category'] == 'INS', 'junction_variant_start'] = var_df['start'] -1
    var_df.loc[var_df['variant_category'] == 'INS', 'junction_variant_stop'] = var_df['start']

    # finish this logic first
    # MNVs
    #var_df.loc[(var_df['variant_category'] == 'MNV') & (var_df['variant_info'].str.contains('{var_df['start']}|{var_df['stop']}')), 'junction_variant_start'] = var_df['start']
    var_df.loc[var_df['variant_category'] == 'MNV', 'junction_variant_stop'] = var_df['stop']

    # format to match junctions
    var_df['variant_info'] = var_df['chromosome_name'] + ':' + var_df['junction_variant_start'].astype('string') + '-' + var_df['junction_variant_stop'].astype('string')

    var_df.to_csv('/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-6-14_indel_fix/H_NJ-HCC1395-HCC1395_annotated.tsv', sep='\t', index=False)
    
    return var_df

#def execute(self):

variant_df = add_junction_coordinates_to_variants('/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-6-14_indel_fix/H_NJ-HCC1395-HCC1395_annotated.tsv') #self.variants
junction_df = pd.read_csv('/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-6-14_indel_fix/H_NJ-HCC1395-HCC1395_filtered.tsv', sep='\t')

merged_df = junction_df.merge(variant_df, on=['variant_info'])

merged_df.to_csv('/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-6-14_indel_fix/H_NJ-HCC1395-HCC1395_combined.tsv', sep='\t', index=False)

# variant_df unique variants = 46

# junction_df unique variants = 13
# merged_df = 12 - missing MNV