from gtfparse import read_gtf
import pandas as pd    

class LoadGtfData():
    def __init__(self, **kwargs):
        self.gtf_file    = kwargs['gtf_file']
        self.output_file = kwargs['output_file']
        self.save_gtf    = kwargs['save_gtf']

    def execute(self):
        print('Converting gtf file to dataframe')

        gtf_df_all = read_gtf(self.gtf_file, usecols=['feature', 'seqname', 'start', 'end', 'transcript_id', 'transcript_biotype', 'transcript_version', 'transcript_support_level', 'exon_number', 'gene_name', 'gene_id'])

        # print('GTF all df')
        # print(gtf_df_all['feature'].unique())
        # print(len(gtf_df_all.index))
        
        # tscript and CDS - this is coding coordinates exon by exon (leaving out 5/3' UTRs in exon body)
        gtf_df = gtf_df_all.loc[(gtf_df_all['feature'].isin(['CDS', 'transcript'])) & (gtf_df_all['transcript_support_level'] == '1') & (gtf_df_all['transcript_biotype'] == 'protein_coding')].replace(["^\s*$"], 'NA', regex=True)
        
        # print('GTF filtered df')
        # print(gtf_df['feature'].unique())
        # print(len(gtf_df.index))
        
        gtf_df = gtf_df.rename(columns={'start': 'cds_start', 'end': 'cds_stop', 'seqname': 'cds_chrom'})
        
        if self.save_gtf:
            gtf_df.to_csv(self.output_file, sep='\t', index=False)

        return gtf_df
