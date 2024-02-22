import pandas as pd
from gtfparse import read_gtf

# can put in docs to run this before pvacsplice
# use to only load gtf_tsv once and not have to include it in all sample runs
# or optional ?

def load_gtf_df(gtf_file):
    print('Converting gtf file to dataframe')

    gtf_df_all = read_gtf(gtf_file, usecols=['feature', 'seqname', 'start', 'end', 'transcript_id', 'transcript_biotype', 'transcript_version', 'transcript_support_level', 'exon_number', 'gene_name', 'gene_id'])

    print('GTF all df')
    print(gtf_df_all['feature'].unique())
    print(len(gtf_df_all.index))
    
    # tscript and CDS - this is coding coordinates exon by exon (leaving out 5/3' UTRs in exon body)
    gtf_df = gtf_df_all.loc[(gtf_df_all['feature'].isin(['CDS', 'transcript'])) & (gtf_df_all['transcript_support_level'] == '1') & (gtf_df_all['transcript_biotype'] == 'protein_coding')].replace(["^\s*$"], 'NA', regex=True)
    
    print('GTF feaure df')
    print(gtf_df['feature'].unique())
    print(len(gtf_df.index))
    
    gtf_df = gtf_df.rename(columns={'start': 'cds_start', 'end': 'cds_stop', 'seqname': 'cds_chrom'})
    
    gtf_df.to_csv(f'{gtf_file}.tsv', sep='\t', index=False)

    # pandas na values = pd.NA
    #gtf_df = pd.read_csv(f'{gtf_file}.tsv', sep='\t')

    #print('Completed')

    #return gtf_df