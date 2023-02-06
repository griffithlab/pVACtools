from gtfparse import read_gtf
import sys

    
# read in gtf file with gtfparse, set index as transcript_id (non-unique index)
gtf_df = read_gtf('/Users/mrichters/Documents/Alt_Splicing/HCC1395/Homo_sapiens.GRCh38.95.gtf', usecols=['feature', 'seqname', 'start', 'end', 'transcript_id', 'transcript_biotype', 'transcript_version', 'transcript_support_level', 'exon_number'])

gtf_df = gtf_df.loc[(gtf_df['feature'] == 'transcript') | (gtf_df['feature'] == 'exon')]

# save file
gtf_df.to_csv('/Users/mrichters/Documents/Alt_Splicing/HCC1395/2022-11-01_gtfparse_test/gtfv95_df.tsv', sep='\t', index=False)

# get size of objects in memory (in bytes)
#print(sys.getsizeof(gtf_df)) # does this actually work? size stays the same even if I slice the data
mem = gtf_df.memory_usage(index=True, deep=True)
total_size = sum(mem)

print(f'{total_size:,} B')
print(f'{total_size/1000000000:.2f} GB')
# using 1.25 GB

# check dtypes of cols (int64: start, end; score: float32; rest: object)
#gtf_df.dtypes

gtf_df_transcripts = gtf_df[(gtf_df['feature'] == 'transcript') & (gtf_df['transcript_id'] == 'ENST00000456328')]
    
# load ensembl dataset
#dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')

# get transcripts info
#df = dataset.query(attributes=["ensembl_transcript_id", "strand", "transcript_start", "transcript_end", "exon_chrom_start", "exon_chrom_end", "genomic_coding_start", "genomic_coding_end"], filters={'link_ensembl_transcript_stable_id': [tscript_id]})
# transcript not found in ensembl (Exception)
#if len(df) == 0:
#    print('Transcript sequence not found in ensembl...skipping')

# drop exons not in coding sequence, sort sequence, reset index w/ coding exons
#df = df.dropna().sort_values(by="Genomic coding start").reset_index(drop=True)            
# convert coordinates from str to int
#df["Genomic coding start"] = df["Genomic coding start"].astype(int)
#df["Genomic coding end"] = df["Genomic coding end"].astype(int)

#return df
