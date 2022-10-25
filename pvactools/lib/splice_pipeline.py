import os
import pandas as pd
import pyfaidx
import shutil
from pvactools.lib.filter_regtools_results import *
from pvactools.lib.junction_to_fasta import *
from pvactools.lib.fasta_to_kmers import *
from pvactools.lib.combine_inputs import *
from pvactools.lib.run_argument_parser import *
from pvactools.lib.input_file_converter import PvacspliceVcfConverter


class JunctionPipeline():
    def __init__(self, **kwargs):
        self.input_file              = kwargs['input_file']
        self.sample_name             = kwargs['sample_name']
        self.output_dir              = kwargs['base_output_dir']
        self.fasta_path              = kwargs['ref_fasta']        
        self.annotated_vcf           = kwargs['annotated_vcf']
        self.input_gtf              = 'insert_cluster_gtf_path'
        self.class_i_epitope_length  = kwargs['class_i_epitope_length']
        self.class_ii_epitope_length = kwargs['class_ii_epitope_length']
        self.class_i_hla             = kwargs['class_i_hla']
        self.class_ii_hla            = kwargs['class_ii_hla']
        self.junction_score          = kwargs.pop('junction_score', 10)
        self.variant_distance        = kwargs.pop('variant_distance', 100)
        self.maximum_transcript_support_level = kwargs.pop('maximum_transcript_support_level', None)
        self.normal_sample_name      = kwargs.pop('normal_sample_name', None)
    
    def execute(self):
        self.filter_regtools_results()
        self.vcf_to_tsv()
        self.create_fastas()
        self.combine_inputs()
        self.junction_to_fasta()
        self.fasta_to_kmers()

    def create_fastas(self):
        fasta_basename = os.path.basename(self.fasta_path)
        alt_fasta_path = os.path.join(self.output_dir, f'{fasta_basename.split(".")[0]}_alt.{".".join(fasta_basename.split(".")[1:])}')
        print("Building alternative fasta")
        size1 = os.path.getsize(self.fasta_path)
        if not os.path.exists(alt_fasta_path):
            shutil.copy(self.fasta_path, alt_fasta_path)
            size2 = os.path.getsize(alt_fasta_path)
            if os.path.exists(alt_fasta_path) and size1 == size2:
                print('Completed')
        elif os.path.exists(alt_fasta_path) and size1 != os.path.getsize(alt_fasta_path):
            print('Fasta transfer is incomplete. Trying again.')
            shutil.copy(self.fasta_path, alt_fasta_path)
            size2 = os.path.getsize(alt_fasta_path)
            if size1 == size2:
                print('Completed')
        elif os.path.exists(alt_fasta_path) and size1 == os.path.getsize(alt_fasta_path):
            print('Alternative fasta already exists. Skipping.')
        print('Creating fasta objects')
        self.personalized_fasta = pyfaidx.FastaVariant(alt_fasta_path, self.annotated_vcf, sample=self.sample_name)
        print('Completed')

    def create_file_path(self, key):
        inputs = {
            'annotated' : '_annotated.tsv',
            'filtered'  : '_filtered.tsv',
            'combined'  : '_combined.tsv',
            'fasta'     : '.transcripts.fa', 
        }
        file_name = os.path.join(self.output_dir, self.sample_name + inputs[key])
        return file_name

    # creates filtered file
    # self.ensembl_version
    def filter_regtools_results(self):
        print('Filtering regtools results')
        filter_params = {
            'input_file'  : self.input_file,
            'output_file' : self.create_file_path('filtered'),
            'score'       : self.junction_score,
            'distance'    : self.variant_distance,
        }
        filter = FilterRegtoolsResults(**filter_params)
        self.filter_df = filter.execute()
        print('Completed')

    # creates annotated file
    def vcf_to_tsv(self):
        convert_params = {
            'input_file'  : self.annotated_vcf,
            'output_file' : self.create_file_path('annotated'),
            'sample_name' : self.sample_name,
        }
        if self.normal_sample_name:
            convert_params['normal_sample_name'] = self.normal_sample_name
        print('Converting .vcf to TSV')
        if os.path.exists(self.create_file_path('annotated')):
            print('Annotated variant file already exists. Skipping.')
        else:
            converter = PvacspliceVcfConverter(**convert_params)
            converter.execute()
            print('Completed')
    
    # creates combined file
    def combine_inputs(self):
        combine_params = {
            'junctions_df'   : self.filter_df,
            'variant_file'   : self.create_file_path('annotated'),
            'sample_name'    : self.sample_name,
            #'output_file'    : self.create_file_path('combined'),
            'maximum_transcript_support_level' : self.maximum_transcript_support_level,
        }
        print('Combining junction and variant information')
        combined = CombineInputs(**combine_params)
        self.combined_df = combined.execute()
        print('Completed')

    # creates transcripts.fa
    def junction_to_fasta(self):
        self.final_combined_df = self.combined_df.copy()
        print('Assembling tumor-specific splicing junctions')   
        for i in self.combined_df.index.unique().to_list():
            print(i)
            junction = self.combined_df.loc[[i], :]         
            for row in junction.itertuples():
                junction_params = {
                    'fasta'          : self.personalized_fasta,
                    'junction_df'    : self.combined_df,
                    'tscript_id'     : row.transcript_name,
                    'chrom'          : row.junction_chrom,
                    'junction_name'  : row.name,
                    'junction_coors' : [row.junction_start, row.junction_stop],
                    'fasta_index'    : row.index,
                    'variant_info'   : row.variant_info,
                    'anchor'         : row.anchor,
                    'strand'         : row.strand,
                    'gene_name'      : row.Gene_name,
                    'output_file'    : self.create_file_path('fasta'),
                    'output_dir'     : self.output_dir,
                    'sample_name'    : self.sample_name,
                    'vcf'            : self.annotated_vcf,
                }
                junctions = JunctionToFasta(**junction_params)
                wt = junctions.create_wt_df()
                if wt.empty:
                    continue
                alt = junctions.create_alt_df()
                if alt.empty:
                    continue
                wt_aa, wt_fs = junctions.get_aa_sequence(wt, 'wt')
                alt_aa, alt_fs = junctions.get_aa_sequence(alt, 'alt')
                if wt_aa == '' or alt_aa == '':
                    print('Amino acid sequence was not produced. Skipping.')
                    continue
                junctions.create_sequence_fasta(wt_aa, alt_aa)
                # df[row, col]
                self.final_combined_df.loc[i, 'wt_protein_length'] = len(wt_aa)
                self.final_combined_df.loc[i, 'alt_protein_length'] = len(alt_aa)
                self.final_combined_df.loc[i, 'frameshift_event'] = alt_fs
        print('Completed')
        self.final_combined_df.to_csv(self.create_file_path('combined'), sep='\t', index=False)
    
    # creates kmer fasta files for input into prediction pipeline
    def fasta_to_kmers(self):
        kmer_params = {
            'fasta'           : self.create_file_path('fasta'),
            'output_dir'      : self.output_dir,
            'class_i_epitope_length' : self.class_i_epitope_length,
            'class_ii_epitope_length': self.class_ii_epitope_length,
            'class_i_hla'     : self.class_i_hla,
            'class_ii_hla'    : self.class_ii_hla,
            'sample_name'     : self.sample_name,
        }           
        fasta = FastaToKmers(**kmer_params)
        fasta.execute()
        print('Completed')
     
        