import os
import sys
import pandas as pd

class FilterRegtoolsResults:
    def __init__(self, **kwargs):
        self.input_file  = kwargs['input_file']
        self.output_file = kwargs['output_file']
        self.gtf_data      = kwargs['gtf_data']
        self.score       = kwargs['score'] #10 reads, 100, 500 -j (cohort stats combined junction tsv filters out junctions with <= 5 reads).  
        self.distance    = kwargs['distance'] #50, 100 bp, 150 -v

    def split_string(self, df, col, delimiter):
        df[col] = df[col].str.split(delimiter)

    def filter_junction_rows(self):
        # open file, rename junction cols for clarity
        junctions = pd.read_csv(self.input_file, sep='\t')
        junctions['transcripts'] = junctions['transcripts'].astype(str)
        junctions = junctions.rename(columns={'chrom':'junction_chrom', 'start':'junction_start', 'end':'junction_stop'})

        # filter on score, strand, and anchor
        filter_junctions = junctions[(junctions['score'] > self.score) & (junctions['strand'] != '?') & (junctions['anchor'].isin(['D', 'A', 'NDA']))].dropna()

        # create variant_start col
        filter_junctions['variant_start'] = filter_junctions['variant_info'].str.split(':|-').str[1].astype('int64')

        # filter by distance: variant_start > start-distance and variant_start < end+distance
        # does strand matter here - no
        final_filter = filter_junctions[
           (filter_junctions['variant_start'] > filter_junctions['junction_start'].astype('int64') - self.distance) & 
           (filter_junctions['variant_start'] < filter_junctions['junction_stop'].astype('int64') + self.distance)
        ].reset_index()

        return final_filter

    def pc_junction_rows(self, filter_junctions):
        # example entry: 0: {'gene_ids': 'ENSG00000122483', 'transcripts': 'ENST00000343253,ENST00000370276,ENST00000401026,ENST00000421014,ENST00000455267'}
        tscript_dict = {i:{'gene_ids': x, 'transcripts': y} for i,(x,y) in enumerate(zip(filter_junctions['gene_ids'], filter_junctions['transcripts']))}

        # filter transcripts by protein_coding and transcript_id
        pc_junctions = pd.DataFrame()

        for k,v in tscript_dict.items():

            # subset df by transcript_id
            gtf_transcripts = self.gtf_data[(self.gtf_data['feature'] == 'transcript') & (self.gtf_data['transcript_id'].isin(v['transcripts'].split(',')))]

            if not gtf_transcripts.empty:
                # add to df
                pc_junctions = pd.concat([pc_junctions, gtf_transcripts])

        # subset of self.gtf_df
        return pc_junctions

    def explode_junction_rows(self, filter_junctions):
        # make transcripts/variants from str to transcript list
        self.split_string(filter_junctions, 'transcripts', ',')
        self.split_string(filter_junctions, 'variant_info', ',')

        # explode the transcript list and variant list
        explode_junctions = filter_junctions.explode('transcripts', ignore_index=True).explode('variant_info', ignore_index=True).drop('index', axis=1) 

        explode_junctions = explode_junctions.rename(columns={'transcripts': 'transcript_id'})

        return explode_junctions

    def merge_and_write(self, pc_junctions, explode_junctions):
        if len(pc_junctions) == 0 & len(explode_junctions) == 0:
            merged_df = pd.DataFrame(columns=[
                'junction_chrom', 'junction_start', 'junction_stop', 'name', 'score', 'strand', 'splice_site',
                'acceptors_skipped', 'exons_skipped', 'donors_skipped', 'anchor', 'known_donor', 'known_acceptor',
                'known_junction', 'transcript_id', 'variant_info', 'feature', 'cds_chrom', 'cds_start', 'cds_stop',
                'transcript_biotype', 'transcript_version', 'transcript_support_level', 'gene_name', 'gene_id'
            ])
        else:
            merged_df = explode_junctions.merge(pc_junctions, on='transcript_id').drop_duplicates()
            merged_df = merged_df.drop(columns=['gene_names', 'gene_ids', 'variant_start', 'exon_number'])
            # drop repetitive or unneeded cols
            # switch strand to numeral
            merged_df['strand'] = merged_df['strand'].replace(['+', '-'], [1, -1])
        # create filtered tsv file
        merged_df.to_csv(self.output_file, sep='\t', index=False, na_rep="NA")

        return merged_df

    def execute(self):
        # filter on score, strand, anchor, and distance; add variant_start column
        filter_junctions = self.filter_junction_rows()

        # filter transcripts by protein_coding and transcript_id
        pc_junctions = self.pc_junction_rows(filter_junctions)

        # explode the transcript list and variant list
        explode_junctions = self.explode_junction_rows(filter_junctions)

        # merge dfs and create associated filtered tsv file
        filtered_df = self.merge_and_write(pc_junctions, explode_junctions)
        return filtered_df
