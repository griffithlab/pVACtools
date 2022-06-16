import os
import pandas as pd
from pybiomart import Dataset

class FilterRegtoolsResults():
    def __init__(self, **kwargs):
        self.input_file  = kwargs['input_file']
        self.output_file = kwargs['output_file']
        self.output_dir  = kwargs['output_dir']
        self.score       = kwargs['score']
        self.distance    = kwargs['distance']

    def execute(self):
        # open file, rename junction cols for clarity
        junctions = pd.read_csv(self.input_file, sep='\t')
        junctions = junctions.rename(columns={'chrom':'junction_chrom', 'start':'junction_start', 'end':'junction_stop'})

        # filter on score, strand, and anchor
        filter_junctions = junctions[(junctions['score'] > self.score) & (junctions['strand'] != '?') & (junctions['anchor'].isin(['D', 'A', 'NDA']))].dropna()
        
        # create variant_start col
        filter_junctions['variant_start'] = filter_junctions['variant_info'].str.split(':|-|,', expand=True)[[1]]
        # filter by distance: variant_start > start-distance and variant_start < end+distance
        # does strand matter here - no
        filter_junctions = filter_junctions[
           (filter_junctions['variant_start'].astype(int) > filter_junctions['junction_start'].astype(int) - self.distance) & 
           (filter_junctions['variant_start'].astype(int) < filter_junctions['junction_stop'].astype(int) + self.distance)
        ].reset_index()

        # example entry: 0: {'gene_ids': 'ENSG00000122483', 'transcripts': 'ENST00000343253,ENST00000370276,ENST00000401026,ENST00000421014,ENST00000455267'}
        tscript_dict = {i:{'gene_ids': x, 'transcripts': y} for i,(x,y) in enumerate(zip(filter_junctions['gene_ids'], filter_junctions['transcripts']))}

        # load ensembl database
        dataset = Dataset(name='hsapiens_gene_ensembl', host='http://useast.ensembl.org')

        # filter transcripts by protein_coding and transcript_id
        pc_junctions = pd.DataFrame()
        # this adds 4 columns to the end of junctions file
        # Transcript_stable_ID, Gene_stable_ID, Gene_name, Transcript_type
        for k,v in tscript_dict.items():
            # return these attributes
            protein_coding = dataset.query(attributes=[
                'external_gene_name',
                'ensembl_gene_id',
                'ensembl_transcript_id', 
                'transcript_biotype'
                ],
                # filter on transcripts and protein_coding
                filters={
                    'link_ensembl_transcript_stable_id': v['transcripts'].split(','),
                    'transcript_biotype': 'protein_coding',
                })
            # add to df
            pc_junctions = pd.concat([pc_junctions, protein_coding])

        # make transcripts from str to transcript list
        filter_junctions['transcripts'] = filter_junctions['transcripts'].str.split(',')
        # explode the transcript list
        explode_junctions = filter_junctions.explode('transcripts', ignore_index=True).drop('index', axis=1)

        # filter - pc_tscripts and unique_all_transcripts - 1 to 1 comparison bc 'transcripts' has been exploded
        pc_junctions = pc_junctions[pc_junctions['Transcript stable ID'].isin(explode_junctions['transcripts'])]

        # merge dfs
        merged_df = explode_junctions.merge(pc_junctions, left_on='transcripts', right_on='Transcript stable ID').drop_duplicates()
        # drop repetitive or unneeded cols
        merged_df = merged_df.drop(columns=['gene_names', 'gene_ids', 'transcripts', 'variant_start', 'Transcript type'])
        # remove spaces from col names
        merged_df.columns = merged_df.columns.str.replace(r'\s+', '_', regex=True)
        # switch strand to numeral
        merged_df['strand'].replace(['+','-'], [1,-1])

        # testing
        # print(len(set(merged_df['name'])))
        # print(len(set(merged_df['Transcript_stable_ID'])))
        # print(len(set(merged_df['Gene_stable_ID'])))

        # create filtered tsv file
        merged_df.to_csv(os.path.join(self.output_dir, self.output_file), sep='\t', index=False)

# debugging
# if __name__ == '__main__':
#     print('Filtering regtools results')
#     filter_params = {
#         'input_file'  : '/Users/mrichters/Documents/Alt_Splicing/HCC1395/junctions.tsv',
#         'output_file' : '/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-6-14_indel_fix/H_NJ-HCC1395-HCC1395_filtered.tsv',
#         'output_dir'  : '/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-6-14_indel_fix',
#         'score'       : 10,
#         'distance'    : 100,
#     }
#     filter = FilterRegtoolsResults(**filter_params)
#     filter.execute()
#     print('Completed')