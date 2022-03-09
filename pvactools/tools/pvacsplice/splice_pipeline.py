import os
import sys
import argparse
from turtle import distance
import pandas as pd
from filter_regtools_results import *
from junction_to_fasta import *
from fasta_to_kmers import *
from pvactools.lib.run_argument_parser import *
from pvactools.lib.input_file_converter import VcfConverter


class JunctionPipeline():
    def __init__(self, **kwargs):
        self.input_file                       = kwargs['input_file']
        self.annotated_vcf                    = kwargs['annotated_vcf']
        self.sample_name                      = kwargs['sample_name']
        self.output_dir                       = kwargs['base_output_dir']
        self.ref_fasta                        = kwargs['ref_fasta']
        self.class_i_epitope_length           = kwargs.pop('class_i_epitope_length', [8,9,10,11])
        self.class_ii_epitope_length          = kwargs.pop('class_ii_epitope_length', [12,13,14,15,16,17,18])
        self.junction_score                   = kwargs.pop('junction_score', 10)
        self.variant_distance                 = kwargs.pop('variant_distance', 100)
        self.maximum_transcript_support_level = kwargs.pop('maximum_transcript_support_level', 1)

    def execute(self):
        self.filter_regtools_results()
        self.vcf_to_tsv()
        self.combine_inputs()
        self.junction_to_fasta()
        self.fasta_to_kmers()

    def create_file_path(self, type):
        inputs = {
            'annotated': '_annotated.tsv',
            'filtered' : '_filtered.tsv',
            'combined' : '_combined.tsv',
            'fasta'    : '_transcripts.fa',
        }
        file_name = self.sample_name + inputs[type]
        return os.path.join(self.output_dir, file_name)

    def filter_regtools_results(self):
        print('Filtering regtools results')
        filter_params = {
            'input_file'  : self.input_file,
            'output_file' : self.create_file_path('filtered'),
            'sample_name' : self.sample_name,
            'output_dir'  : self.output_dir,
            'score'       : self.junction_score,
            'distance'    : self.variant_distance,
            'tsl'         : self.maximum_transcript_support_level,
        }
        filter = FilterRegtoolsResults(**filter_params)
        filter.execute()
        print('Completed')

    def vcf_to_tsv(self):
        print('Converting .vcf to TSV')
        convert_params = {
            'input_file'  : self.annotated_vcf,
            'output_file' : self.create_file_path('annotated'),
            'sample_name' : self.sample_name,
        }
        converter = VcfConverter(**convert_params)
        converter.execute()
        print('Completed')
    
    def combine_inputs(self):
        print('Combining junction and variant information')
        # should always have the same column names (specify regtools version)
        junctions = pd.read_csv(self.create_file_path('filtered'), sep='\t')
        # dif vcfs might have dif column names? (specify in docs)
        annotated = pd.read_csv(self.create_file_path('annotated'), sep='\t')
        annotated.rename(columns={
            'chromosome_name' : 'variant_chrom', 
            'start'           : 'variant_start', 
            'stop'            : 'variant_stop',
            }, inplace=True
        )
        # remove version number to compare with filtered junctions file
        annotated['transcript_name'] = annotated['transcript_name'].str.split('.', expand=True)[[0]]
        # merge dfs by variant_info, gene_id, gene_name, transcript_id
        self.merged_df = junctions.merge(annotated, on=[
            'variant_chrom', 
            'variant_start', 
            'variant_stop', 
            'transcript_name', 
            'ensembl_gene_id', 
            'gene_name',
            ]
        )
        # uniquify tsv indexes
        self.merged_df['index'] = self.merged_df['index'] + '.' + self.merged_df['name'] 
        self.merged_df.to_csv(self.create_file_path('combined'), sep='\t', index=False)
        print('Completed')


    def junction_to_fasta(self):
        print('Assembling tumor-specific splicing junctions')
        filtered_df = pd.read_csv(self.create_file_path('combined'), sep='\t')
        for i in filtered_df.index.unique().to_list():
            junction = filtered_df.loc[[i], :]
            for row in junction.itertuples():
                junction_params = {
                    'fasta_path'     : self.ref_fasta,
                    'tscript_id'     : row.transcript_name,
                    'chrom'          : row.junction_chrom,
                    'junction_name'  : row.name,
                    'junction_coors' : [row.junction_start, row.junction_stop],
                    'index'          : row.index, 
                    'anchor'         : row.anchor,
                    'strand'         : row.strand,
                    'gene_name'      : row.gene_name,
                    'output_file'    : self.create_file_path('fasta'),
                    'output_dir'     : self.output_dir,
                }
                junctions = JunctionToFasta(**junction_params)
                wt = junctions.create_wt_df()
                if wt.empty:
                    continue
                mut = junctions.create_alt_df()
                if mut.empty:
                    continue
                wt_aa = junctions.get_aa_sequence(wt)
                if not wt_aa:
                    continue
                mut_aa = junctions.get_aa_sequence(mut)
                junctions.create_sequence_fasta(wt_aa, mut_aa)
        print('Completed')
    

    def fasta_to_kmers(self):
        print('Creating a fasta of kmer peptides from altered portion of transcripts')
        kmer_params = {
            'fasta'           : self.create_file_path('fasta'),
            'output_dir'      : self.output_dir,
            'epitope_lengths' : self.class_i_epitope_length + self.class_ii_epitope_length, 
            'combined_df'    : self.merged_df,
        }
        fasta = FastaToKmers(**kmer_params)
        fasta.execute()
        print('Completed')

