import gtfparse
import sys


class LoadGtfData:
    def __init__(self, **kwargs):
        self.gtf_file    = kwargs['gtf_file']
        self.output_file = kwargs['output_file']
        self.save_gtf    = kwargs['save_gtf']  # default false
        self.biotypes    = kwargs.pop('biotypes', ['protein_coding'])

    def execute(self):

        # make sure running pandas not polars
        gtf_df_all = gtfparse.read_gtf(self.gtf_file, usecols=[
            'feature', 'seqname', 'start', 'end', 'transcript_id',
            'transcript_biotype', 'transcript_version', 'transcript_support_level', 'exon_number', 'gene_name',
            'gene_id'], result_type='pandas')

        if 'transcript_support_level' in gtf_df_all:
            gtf_df_all['transcript_support_level'] = gtf_df_all['transcript_support_level'].str.replace(r' (assigned to previous version \d)', '', regex=True)
            gtf_df_all['transcript_support_level'] = gtf_df_all['transcript_support_level'].str.replace(r'^$', 'NA', regex=True)
        else:
            gtf_df_all['transcript_support_level'] = 'Not Supported'

        # transcript and CDS - this is coding coordinates exon by exon (leaving out 5/3' UTRs in exon body)
        if 'transcript_biotype' in gtf_df_all:
            gtf_df = gtf_df_all.loc[
                (gtf_df_all['feature'].isin(['CDS', 'transcript'])) &
                (gtf_df_all['transcript_biotype'].isin(self.biotypes))
                ].replace(["^\s*$"], 'NA', regex=True)
        else:
            gtf_df = gtf_df_all.loc[
                (gtf_df_all['feature'].isin(['CDS', 'transcript']))
                ].replace(["^\s*$"], 'NA', regex=True)
            gtf_df['transcript_biotype'] = 'Not Supported'

        gtf_df = gtf_df.rename(columns={'start': 'cds_start', 'end': 'cds_stop', 'seqname': 'cds_chrom'})

        if gtf_df.empty:
            sys.exit('The GTF dataframe is empty. Please check your input file.')

        if self.save_gtf:
            gtf_df.to_csv(self.output_file, sep='\t', index=False)
            print(f'GTF TSV output file has been saved as: {self.output_file}')

        return gtf_df
