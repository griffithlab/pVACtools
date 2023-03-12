from pvactools.lib.filter_regtools_results import *
from pvactools.lib.junction_to_fasta import *
from pvactools.lib.fasta_to_kmers import *
from pvactools.lib.combine_inputs import *
from pvactools.lib.run_argument_parser import *
from pvactools.lib.input_file_converter import PvacspliceVcfConverter
from pvactools.lib.load_gtf_data import *


class JunctionPipeline:
    def __init__(self, **kwargs):
        self.input_file              = kwargs['input_file']
        self.sample_name             = kwargs['sample_name']
        self.output_dir              = kwargs['base_output_dir']
        self.fasta_path              = kwargs['ref_fasta']        
        self.annotated_vcf           = kwargs['annotated_vcf']
        self.gtf_file                = kwargs['gtf_file']
        self.class_i_epitope_length  = kwargs['class_i_epitope_length']
        self.class_ii_epitope_length = kwargs['class_ii_epitope_length']
        self.class_i_hla             = kwargs['class_i_hla']
        self.class_ii_hla            = kwargs['class_ii_hla']
        self.junction_score          = kwargs['junction_score']
        self.variant_distance        = kwargs['variant_distance']
        self.tsl                     = kwargs['maximum_transcript_support_level']
        self.normal_sample_name      = kwargs.pop('normal_sample_name', None)
        self.save_gtf                = kwargs['save_gtf']

    @staticmethod
    def file_exists(file_path: str, file_type: str):
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            print(f"{file_type} file already exists. Skipping.")
            exists = True
        else:
            exists = False
        return exists

    def execute(self):
        self.load_gtf_data()
        self.vcf_to_tsv()
        self.junction_to_fasta()
        self.fasta_to_kmers()

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
            'gtf_file'    : self.gtf_file,
            'output_file' : self.create_file_path('gtf'),
            'save_gtf'    : self.save_gtf, # default no but option to save for running cohorts processed with the same reference data
            'tsl'         : self.tsl,
            }
            gtf_data = LoadGtfData(**gtf_params)
            gtf_df = gtf_data.execute()
            print('Completed')
        return gtf_df

    def create_file_path(self, key):
        inputs = {
            'gtf'       : '_gtf.tsv',
            'annotated' : '_annotated.tsv',
            'filtered'  : '_filtered.tsv',
            'combined'  : '_combined.tsv',
            'fasta'     : '.transcripts.fa',
        }
        file_name = os.path.join(self.output_dir, self.sample_name + inputs[key])
        return file_name

    def filter_regtools_results(self):
        if self.file_exists(self.create_file_path('filtered'), 'Filtered'):
            filter_df = pd.read_csv(self.create_file_path('filtered'), sep='\t')
        else:
            print('Filtering regtools results')
            filter_params = {
                'input_file'  : self.input_file,
                'output_file' : self.create_file_path('filtered'),
                'gtf_data'    : self.gtf_df,
                'score'       : self.junction_score,
                'distance'    : self.variant_distance,
            }
            filter_object = FilterRegtoolsResults(**filter_params)
            filter_df = filter_object.execute()
            print('Completed')
        return filter_df

    def vcf_to_tsv(self):
        if self.file_exists(self.create_file_path('annotated'), 'VCF TSV'):
            pass
        else:
            print('Converting VCF to TSV')
            convert_params = {
                'input_file'  : self.annotated_vcf,
                'output_file' : self.create_file_path('annotated'),
                'sample_name' : self.sample_name,
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
            combine_params = {
                'junctions_df' : self.filter_regtools_results(),
                'variant_file' : self.create_file_path('annotated'),
                'output_dir'   : self.output_dir,
                'output_file'  : self.create_file_path('combined'),
            }
            print('Merging junction and variant info')
            combined = CombineInputs(**combine_params)
            combined_df = combined.execute()
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
                print(i)
                junction = combined_df.loc[[i], :]
                for row in junction.itertuples():
                    junction_params = {
                        'fasta_path'     : self.fasta_path,
                        'junction_df'    : combined_df,
                        'gtf_df'         : self.load_gtf_data(),
                        'tscript_id'     : row.transcript_id,
                        'chrom'          : row.junction_chrom,
                        'junction_name'  : row.name,
                        'junction_coors' : [row.junction_start, row.junction_stop],
                        'fasta_index'    : row.index,
                        'variant_info'   : row.variant_info,
                        'anchor'         : row.anchor,
                        'strand'         : row.strand,
                        'gene_name'      : row.gene_name,
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
            # update file
            combined_df.to_csv(self.create_file_path('combined'), sep='\t', index=False)
            print('Completed')

    def fasta_to_kmers(self):
        if self.file_exists(self.create_file_path('fasta'), 'Peptide fasta'):
            pass
        else:
            print('Generating peptides from novel junction sequences')
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
     
        