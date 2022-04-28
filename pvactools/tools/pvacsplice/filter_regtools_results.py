import os
import sys
import pandas as pd
from pybiomart import Dataset
import argparse

class FilterRegtoolsResults():
    def __init__(self, **kwargs):
        self.input_file  = kwargs['input_file']
        self.output_file = kwargs['output_file']
        #self.sample_name = kwargs['sample_name']
        self.output_dir  = kwargs['output_dir']
        self.score       = kwargs['score']
        self.distance    = kwargs['distance']
        self.tsl         = kwargs['tsl']

    def execute(self):
        # open file, rename junction cols for clarity
        junctions = pd.read_csv(self.input_file, sep='\t')
        junctions = junctions.rename(columns={'chrom':'junction_chrom', 'start':'junction_start', 'end':'junction_stop'})

        # filter on score, strand, and anchor
        filter_junctions = junctions[(junctions['score'] > self.score) & (junctions['strand'] != '?') & (junctions['anchor'].isin(['D', 'A', 'NDA']))].dropna().reset_index()
        # split variant_info into sep. cols
        filter_junctions[['variant_chrom', 'variant_start', 'variant_stop']] = filter_junctions['variant_info'].str.split(':|-|,', expand=True)[[0, 1, 2]]

        # filter by distance: variant_start > start-distance and variant_start < end+distance
        # does strand matter here - no
        filter_junctions = filter_junctions[
           (filter_junctions['variant_start'].astype(int) > filter_junctions['junction_start'].astype(int) - self.distance) & 
           (filter_junctions['variant_start'].astype(int) < filter_junctions['junction_stop'].astype(int) + self.distance)
        ] 
        
        # example entry: 0: {'gene_ids': 'ENSG00000122483', 'transcripts': 'ENST00000343253,ENST00000370276,ENST00000401026,ENST00000421014,ENST00000455267'}
        tscript_dict = {i:{'gene_ids': x, 'transcripts': y} for i,(x,y) in enumerate(zip(filter_junctions['gene_ids'], filter_junctions['transcripts']))}

        # load ensembl database
        dataset = Dataset(name='hsapiens_gene_ensembl', host='http://useast.ensembl.org')

        # filter transcripts by protein_coding and there is a tsl assigned to tscript
        pc_junctions = pd.DataFrame()
        for k,v in tscript_dict.items():
            # return these attributes
            protein_coding = dataset.query(attributes=[
                'ensembl_transcript_id', 
                'ensembl_gene_id', 
                'external_gene_name',
                'transcript_tsl',
                'transcript_biotype'
                ],
                # filter on transcripts in tscript_dict
                filters={
                    'link_ensembl_transcript_stable_id': v['transcripts'].split(','),
                    'transcript_biotype': 'protein_coding',
                    'transcript_tsl': True,
                })
            if self.tsl:
                # filter by tsl=1
                protein_coding = protein_coding[protein_coding['Transcript support level (TSL)'].str.contains(f'tsl1')]
            # add to df
            pc_junctions = pc_junctions.append(protein_coding)

        # save file
        #pc_junctions.to_csv('/Users/mrichters/Desktop/Alt_Splicing/HCC1395/2022-3-29_pvacsplice/transcript_filtering.tsv', sep='\t')

        # now choose junctions based on remaining transcripts:

        # make transcripts from str to list and explode the list
        filter_junctions['transcripts'] = filter_junctions['transcripts'].str.split(',')
        explode_junctions = filter_junctions.explode('transcripts', ignore_index=True).drop('index', axis=1)

        #explode_junctions.to_csv('/Users/mrichters/Desktop/Alt_Splicing/HCC1395/2022-3-29_pvacsplice/explode_transcripts.tsv', sep='\t')

        # filter explode_junctions with pc['transcripts'] list
        final_df = explode_junctions[explode_junctions['transcripts'].isin(pc_junctions['Transcript stable ID'])]
        # merge and drop duplicates
        final_df = final_df.merge(pc_junctions, left_on='transcripts', right_on='Transcript stable ID').drop_duplicates()
        # remove spaces from col names and switch strand to numeral
        final_df.columns = final_df.columns.str.replace(r'\s+', '_', regex=True)
        final_df['strand'] = final_df['strand'].replace(['+','-'], [1,-1])

        final_df.to_csv(os.path.join(self.output_dir, self.output_file), sep='\t', index=False)