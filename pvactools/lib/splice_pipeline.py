import shutil
import os
import pandas as pd
from pvactools.lib.filter_regtools_results import FilterRegtoolsResults
from pvactools.lib.junction_to_fasta import JunctionToFasta
from pvactools.lib.fasta_to_kmers import FastaToKmers
from pvactools.lib.combine_inputs import CombineInputs
from pvactools.lib.input_file_converter import PvacspliceVcfConverter
from pvactools.lib.load_gtf_data import LoadGtfData

class JunctionPipeline:
    def __init__(self, **kwargs):
        self.input_file = kwargs['input_file']
        self.sample_name = kwargs['sample_name']
        # default pvacsplice output dir (where splice files are currently)
        self.output_dir = kwargs['junctions_dir']
        self.fasta_path = kwargs['ref_fasta']
        self.annotated_vcf = kwargs['annotated_vcf']
        self.pass_only = kwargs['pass_only']
        self.gtf_file = kwargs['gtf_file']
        self.class_i_epitope_length = kwargs['class_i_epitope_length']
        self.class_ii_epitope_length = kwargs['class_ii_epitope_length']
        self.class_i_hla = kwargs['class_i_hla']
        self.class_ii_hla = kwargs['class_ii_hla']
        self.junction_score = kwargs['junction_score']
        self.variant_distance = kwargs['variant_distance']
        self.normal_sample_name = kwargs.pop('normal_sample_name', None)
        self.save_gtf = kwargs['save_gtf']
        self.keep_tmp_files = kwargs['keep_tmp_files']
        self.biotypes = kwargs['biotypes']
        self.gtf_data = self.load_gtf_data()
        self.tmp_dir = os.path.join(self.output_dir, 'tmp')
        os.makedirs(self.tmp_dir, exist_ok=True)

    @staticmethod
    def file_exists(file_path: str, file_type: str):
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            print(f"{file_type} file already exists. Skipping.")
            exists = True
        else:
            exists = False
        return exists

    def execute(self):
        self.vcf_to_tsv()
        self.junction_to_fasta()
        self.fasta_to_kmers()

# if self.save_gtf --> should not have a gtf tsv file present
# if self.save_gtf == True:
    # save gtf file in output dir
# elif False:
    # write gtf_file to tmp_dir

    def load_gtf_data(self):
        # option 1: load from preexisting gtf tsv
        # option 2: create gtf_df and DON'T save automatically
        # option 3: create gtf_df and DO save automatically
        if self.file_exists(self.create_file_path('gtf'), 'GTF TSV'):
            # load from tsv file
            gtf_df = pd.read_csv(self.create_file_path('gtf'), sep='\t')
        else:
            print('Converting GTF file to TSV')
            gtf_params = {
                'gtf_file': self.gtf_file,
                'output_file': self.create_file_path('gtf'),
                'save_gtf': self.save_gtf,
                # default no but option to save for running cohorts processed with the same reference data
                'biotypes': self.biotypes,
            }
            gtf_data = LoadGtfData(**gtf_params)
            gtf_df = gtf_data.execute()
            print('Completed')
        return gtf_df

    def create_file_path(self, key, temp=False):
        inputs = {
            'gtf': '_gtf.tsv',
            'annotated': '_annotated.tsv',
            'filtered': '_filtered.tsv',
            'combined': '_combined.tsv',
            'fasta': '.transcripts.fa',
        }
        if temp:
            file_name = os.path.join(self.tmp_dir, self.sample_name + inputs[key])
        else:
            file_name = os.path.join(self.output_dir, self.sample_name + inputs[key])

        return file_name

    def filter_regtools_results(self):
        if self.file_exists(self.create_file_path('filtered', temp=True), 'Filtered'):
            filter_df = pd.read_csv(self.create_file_path('filtered', temp=True), sep='\t')
        else:
            print('Filtering regtools results')
            filter_params = {
                'input_file': self.input_file,
                'output_file': self.create_file_path('filtered', temp=True),
                'gtf_data': self.gtf_data,
                'score': self.junction_score,
                'distance': self.variant_distance,
            }
            filter_object = FilterRegtoolsResults(**filter_params)
            filter_df = filter_object.execute()
            print('Completed')
        return filter_df

    def vcf_to_tsv(self):
        if self.file_exists(self.create_file_path('annotated', temp=True), 'VCF TSV'):
            pass
        else:
            print('Converting VCF to TSV')
            convert_params = {
                'input_file': self.annotated_vcf,
                'pass_only': self.pass_only,
                'output_file': self.create_file_path('annotated', temp=True),
                'sample_name': self.sample_name,
            }
            if self.normal_sample_name:
                convert_params['normal_sample_name'] = self.normal_sample_name
            converter = PvacspliceVcfConverter(**convert_params)
            converter.execute()
            print('Completed')

    def combine_inputs(self):
        if self.file_exists(self.create_file_path('combined'), 'Combined junction'):
            combined_df = pd.read_csv(self.create_file_path('combined'), sep='\t')
        else:
            print('Merging junction and variant info')
            junctions_df = self.filter_regtools_results()
            if len(junctions_df) == 0:
                sys.exit("The RegTools junctions TSV file doesn't contain any splice sites supported by pVACsplice. Aborting.")

            combine_params = {
                'junctions_df': junctions_df,
                'variant_file': self.create_file_path('annotated', temp=True),
                'output_dir': self.output_dir,
                'output_file': self.create_file_path('combined'),
            }
            combined = CombineInputs(**combine_params)
            combined_df = combined.execute()
            if len(combined_df) == 0:
                sys.exit("Combined dataset is empty. Aborting.")
            print('Completed')
        return combined_df

    # creates transcripts.fa
    def junction_to_fasta(self):
        combined_df = self.combine_inputs()
        if self.file_exists(self.create_file_path('fasta'), 'Junction fasta'):
            pass
        else:
            print('Assembling tumor-specific splicing junctions')
            for i in combined_df.index.to_list():
                junction = combined_df.loc[[i], :]
                for row in junction.itertuples():
                    junction_params = {
                        'fasta_path': self.fasta_path,
                        'junction_df': combined_df,
                        'gtf_df': self.gtf_data,
                        'tscript_id': row.transcript_id,
                        'chrom': row.junction_chrom,
                        'junction_name': row.name,
                        'junction_coors': [row.junction_start, row.junction_stop],
                        'fasta_index': row.index,
                        'variant_info': row.variant_info,
                        'anchor': row.anchor,
                        'strand': row.strand,
                        'gene_name': row.gene_name,
                        'output_file': self.create_file_path('fasta'),
                        'output_dir': self.output_dir,
                        'sample_name': self.sample_name,
                        'vcf': self.annotated_vcf,
                    }
                    junctions = JunctionToFasta(**junction_params)
                    wt = junctions.create_wt_df()
                    if wt.empty:
                        continue
                    alt = junctions.create_alt_df()
                    if alt.empty:
                        continue
                    wt_aa, wt_fs = junctions.get_aa_sequence(wt)
                    alt_aa, alt_fs = junctions.get_aa_sequence(alt)
                    if wt_aa == '' or alt_aa == '':
                        print('No amino acid sequence was produced. Skipping.')
                        continue
                    # creates output transcript fasta
                    junctions.create_sequence_fasta(wt_aa, alt_aa)
                    # df[row, col]
                    combined_df.loc[i, 'wt_protein_length'] = len(wt_aa)
                    combined_df.loc[i, 'alt_protein_length'] = len(alt_aa)
                    combined_df.loc[i, 'frameshift_event'] = alt_fs
            combined_df = combined_df.dropna(subset=['wt_protein_length', 'alt_protein_length', 'frameshift_event'])
            if len(combined_df) == 0:
                sys.exit("Unable to determine junction peptide sequence for any of the candidates. Aborting.")
            # update file
            combined_df.to_csv(self.create_file_path('combined'), sep='\t', index=False, na_rep="NA")
            print('Completed')

    def fasta_to_kmers(self):
        for el in self.class_i_epitope_length + self.class_ii_epitope_length:
            fasta_file = f'{self.output_dir}.{self.sample_name}.{el}.fa'
            if os.path.exists(fasta_file):
                print(f'{el}mer fasta already exists. Skipping.')
                continue
            else:
                print(f'Generating {el}mer peptides from novel junction sequences')
                kmer_params = {
                    'fasta': self.create_file_path('fasta'),
                    'output_dir': self.tmp_dir,
                    'class_i_epitope_length': self.class_i_epitope_length,
                    'class_ii_epitope_length': self.class_ii_epitope_length,
                    'class_i_hla': self.class_i_hla,
                    'class_ii_hla': self.class_ii_hla,
                    'sample_name': self.sample_name,
                }
                fasta = FastaToKmers(**kmer_params)
                fasta.execute()
                print('Completed')
