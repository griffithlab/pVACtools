import os
import sys
import pandas as pd
from pybiomart import Dataset
import argparse

class FilterRegtoolsResults():
    def __init__(self, **kwargs):
        self.input_file  = kwargs['input_file']
        self.output_file = kwargs['output_file']
        self.sample_name = kwargs['sample_name']
        self.output_dir  = kwargs['output_dir']
        self.score       = kwargs['score']
        self.distance    = kwargs['distance']
        self.tsl         = kwargs['tsl']


    def execute(self):
        # load in ensembl dataset
        dataset = Dataset(name='hsapiens_gene_ensembl', host='http://useast.ensembl.org')
        
        # read in regtools junctions.tsv
        junctions = pd.read_csv(self.input_file, sep='\t')
        junctions = junctions.rename(columns={'chrom':'junction_chrom', 'start':'junction_start', 'end':'junction_stop'})

        # filter junctions by score, strand, anchor, remove NAs, reset index
        filter_junctions = junctions[(junctions['score'] > self.score) & (junctions['strand'] != '?') & (junctions['anchor'].isin(['D', 'A', 'NDA']))].dropna().reset_index()
        if len(filter_junctions) == 0:
            raise ValueError('No junctions passed score and anchor filters. Exiting...')
            
        # expand variant_info field into chrom, start, end
        filter_junctions[['variant_chrom', 'variant_start', 'variant_stop']] = filter_junctions['variant_info'].str.split(':|-|,', expand=True)[[0, 1, 2]]

        # filter by chrom (necessary ??)
        filter_junctions[filter_junctions['junction_chrom'] == filter_junctions['variant_chrom']]

        # filter by distance: variant_start > start-distance and variant_start < end+distance 
        # double check this logic
        filter_junctions = filter_junctions[
            (filter_junctions['variant_start'].astype(int) > filter_junctions['junction_start'].astype(int) - self.distance) & 
            (filter_junctions['variant_start'].astype(int) < filter_junctions['junction_stop'].astype(int) + self.distance)
        ] 

        # get transcript info
        tscript_dict = {i:{'gene_ids': x, 'transcripts': y} for i,(x,y) in enumerate(zip(filter_junctions['gene_ids'], filter_junctions['transcripts']))}

        # filter by transcript biotype
        pc_junctions = pd.DataFrame()
        for k,v in tscript_dict.items():
            protein_coding = dataset.query(attributes=[
                'ensembl_transcript_id', 
                'ensembl_gene_id', 
                'external_gene_name',
                'transcript_tsl',
                ], 
                filters={
                    'link_ensembl_transcript_stable_id': v['transcripts'].split(','),
                    'transcript_biotype': 'protein_coding',
                    'transcript_tsl': True,
                })
            if self.tsl:
                protein_coding = protein_coding[protein_coding['Transcript support level (TSL)'].str.contains(f'tsl1')]
            pc_junctions = pc_junctions.append(protein_coding)
            tscript_dict[k]['transcript_name'] = ','.join(protein_coding['Transcript stable ID'].unique())
            tscript_dict[k]['ensembl_gene_id'] = ','.join(protein_coding['Gene stable ID'].unique())
            tscript_dict[k]['gene_name'] = ','.join(protein_coding['Gene name'].unique())

        # merge new values with filter_junctions
        tscript_df = pd.DataFrame(tscript_dict).transpose().replace('', float('NaN')).dropna().drop_duplicates()
        merged_junctions = filter_junctions.merge(tscript_df, on=['transcripts', 'gene_ids'])

        # filter merged_junctions to include fields for junction_to_fasta.py
        merged_junctions['strand'] = merged_junctions['strand'].replace(['+','-'], [1,-1])
        final_merged = merged_junctions.drop('index', axis=1) #[['name', 'pc_transcripts','chrom','start','end','anchor','strand','pc_gene_names']]
        
        # expand tscripts to one per line
        pd.options.mode.chained_assignment = None
        final_merged['transcript_name'] = final_merged.transcript_name.apply(lambda x: x.split(','))
        
        final_merged = final_merged.explode('transcript_name')

        final_merged.to_csv(os.path.join(self.output_dir, self.output_file), index=False, sep='\t')

        return final_merged
