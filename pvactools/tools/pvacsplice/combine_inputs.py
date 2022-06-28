import os
import pandas as pd

class CombineInputs():
    def __init__(self, **kwargs):
        self.junctions   = kwargs['junctions_file']
        self.variants    = kwargs['variant_file']
        self.sample_name = kwargs['sample_name']
        self.output_dir  = kwargs['output_dir']
        self.tsl         = kwargs['tsl']
        self.output_prefix = f'{self.output_dir}/{self.sample_name}'

    def add_junction_coordinates_to_variants(self):
        # read in df
        var_df = pd.read_csv(self.variants, sep='\t')

        # filter variants by tsl and biotype
        var_df = var_df[(var_df['transcript_support_level'] <= self.tsl) & (var_df['biotype'] == 'protein_coding')]

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

        var_df.to_csv(f'{self.output_prefix}_annotated_filtered.tsv', sep='\t', index=False)
        
        return var_df

    def execute(self):
        variant_df = self.add_junction_coordinates_to_variants()
        # print(set(variant_df['variant_info']))

        junction_df = pd.read_csv(f'{self.output_prefix}_filtered.tsv', sep='\t')
        # print(set(junction_df['variant_info']))

        # merge by transcript and variant coors
        merged_df = junction_df.merge(variant_df, on=['transcript_name', 'variant_info']).drop_duplicates()

        # create index to match with kmers
        merged_df['fasta_index'] = merged_df['Gene_name'] + '.' + merged_df['transcript_name'] + '.' + merged_df['name'] + '.' + merged_df['variant_info'] + '.' + merged_df['anchor']
        
        merged_df.to_csv(f'{self.output_prefix}_combined.tsv', sep='\t', index=False)
        # print(set(merged_df['variant_info']))

# debugging
if __name__ == '__main__':
    print('Combine junction and variant information')
    combine_params = {
        'junctions_file' : '/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-6-20_first_test/H_NJ-HCC1395-HCC1395_filtered.tsv',
        'variant_file'   : '/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-6-20_first_test/H_NJ-HCC1395-HCC1395_annotated.tsv',
        'sample_name'    : 'H_NJ-HCC1395-HCC1395',
        'output_dir'     : '/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-6-20_first_test',
        'tsl' : 1,
    }
    combined = CombineInputs(**combine_params)
    combined.execute()
    print('Completed')