import os
import sys
import argparse
import pandas as pd
from filter_regtools_results import *
from junction_to_fasta import *
from fasta_to_kmers import *
from lib.run_argument_parser import *


class JunctionPipeline():
    def __init__(self, **kwargs):
        self.input_file                       = kwargs['input_file']
        self.sample_name                      = kwargs['sample_name']
        self.base_output_dir                  = kwargs['base_output_dir']
        self.ref_fasta                        = kwargs['ref_fasta']
        self.class_i_epitope_length           = kwargs['class_i_epitope_length']
        self.class_ii_epitope_length          = kwargs['class_ii_epitope_length']
        self.junction_score                   = kwargs['junction_score']
        self.variant_distance                 = kwargs['variant_distance']
        self.maximum_transcript_support_level = kwargs['maximum_transcript_support_level']

    def execute(self):
        self.filter_regtools_results()
        self.junction_to_fasta()
        self.convert_fasta_to_peptides()

    def filter_regtools_results(self):
        print('Filtering regtools results')
        FilterRegtoolsResults(
            input_file=self.input_file,
            sample_name=self.sample_name, 
            output_dir=self.base_output_dir, 
            score=self.junction_score,
            distance=self.variant_distance, 
            tsl=self.maximum_transcript_support_level
            ).filter_regtools_results()
        print('Completed')

    def junction_to_fasta(self):
        print('Converting tumor-specific splicing events to alternate transcript')
        filtered_df = pd.read_csv(f'{self.base_output_dir}/{self.sample_name}_filtered_junctions.tsv', sep='\t')
        for i in filtered_df.index.unique().to_list():
            junction = filtered_df.loc[[i], :]
            for row in junction.itertuples():
                j = JunctionToFasta(
                    fasta_path=self.ref_fasta,  
                    tscript_id=row.pc_transcripts, 
                    chrom=row.chrom, 
                    junction_name=row.name,
                    junction_coors=[row.start, row.end], 
                    anchor=row.anchor, 
                    strand=row.strand, 
                    gene_name=row.pc_gene_names,
                    output_dir=self.base_output_dir
                    )
                wt = j.createWtDataframe()
                if wt.empty:
                    continue
                mut = j.createMtDataframe()
                if mut.empty:
                    continue
                wt_aa = j.findProtein(wt)
                if not wt_aa:
                    continue
                mut_aa = j.findProtein(mut)
                j.createFasta(wt_aa, mut_aa)
        print('Completed')
    

    def convert_fasta_to_peptides(self):
        print('Creating a fasta of kmer peptides from altered portion of transcripts')
        p = FastaToKmers(
            tscript_fasta=f'{self.base_output_dir}/assembled_transcripts.fa',
            output_dir=self.base_output_dir,
            classI_lengths=self.class_i_epitope_length,
            classII_lengths=self.class_ii_epitope_length,
            )
        p.loop_through_tscripts()
        p.create_fasta_pvacbind()
        print('Completed')

