from pybiomart import Dataset

def load_ensembl_data(tscript_id:str):
    
    # load ensembl dataset
    dataset = Dataset(name='hsapiens_gene_ensembl', host='http://useast.ensembl.org')
    
    # get transcripts info
    df = dataset.query(attributes=["ensembl_transcript_id", "strand", "transcript_start", "transcript_end", "exon_chrom_start", "exon_chrom_end", "genomic_coding_start", "genomic_coding_end"], filters={'link_ensembl_transcript_stable_id': [tscript_id]})
    # transcript not found in ensembl (Exception)
    if len(df) == 0:
        print('Transcript sequence not found in ensembl...skipping')
    
    # drop exons not in coding sequence, sort sequence, reset index w/ coding exons
    df = df.dropna().sort_values(by="Genomic coding start").reset_index(drop=True)            
    # convert coordinates from str to int
    df["Genomic coding start"] = df["Genomic coding start"].astype(int)
    df["Genomic coding end"] = df["Genomic coding end"].astype(int)

    return df
