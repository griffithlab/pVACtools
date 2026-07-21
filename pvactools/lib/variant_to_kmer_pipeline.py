import os

from pvactools.lib.input_to_kmer_pipeline import InputToKmerPipeline
from pvactools.lib.variant_to_fasta import VariantToFasta
from pvactools.lib.fasta_to_kmers import VariantFastaToKmers
from pvactools.lib.input_file_converter import VcfConverter

class VariantToKmerPipeline(InputToKmerPipeline):
    def __init__(self, **kwargs):
        self.input_file = kwargs['input_file']
        self.output_dir = kwargs['output_dir']
        self.pass_only = kwargs.pop('pass_only', False)
        self.sample_name = kwargs.pop('sample_name', 'tmp')
        self.normal_sample_name = kwargs.pop('normal_sample_name', None)
        self.proximal_variants_vcf = kwargs.pop('proximal_variants_vcf', None)
        self.biotypes = kwargs.pop('biotypes', ['protein_coding'])
        self.allow_incomplete_transcripts = kwargs.pop('allow_incomplete_transcripts', False)
        self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', 1000)
        self.class_i_epitope_length = kwargs.pop('class_i_epitope_length', None)
        self.class_ii_epitope_length = kwargs.pop('class_ii_epitope_length', None)
        self.class_i_hla = kwargs.pop('class_i_hla', None)
        self.class_ii_hla = kwargs.pop('class_ii_hla', None)
        self.flanking_bases = kwargs.pop('flanking_bases', None)
        if self.flanking_bases is None:
            self.flanking_bases = max(self.choose_final_lengths())

    def create_file_path(self, key):
        inputs = {
            'tsv': '.tsv',
            'proximal_variants_tsv': '.proximal_variants.tsv',
            'fasta': '.transcripts.fa',
        }
        file_name = os.path.join(self.output_dir, self.sample_name + inputs[key])

        return file_name

    def input_to_tsv(self):
        if self.file_exists(self.create_file_path('tsv'), 'TSV'):
            pass
        else:
            print("Converting VCF to TSV")
            params = {
                'input_file' : self.input_file,
                'output_file': self.create_file_path('tsv'),
                'pass_only'  : self.pass_only,
                'sample_name': self.sample_name,
                'normal_sample_name': self.normal_sample_name,
                'proximal_variants_vcf': self.proximal_variants_vcf,
                'proximal_variants_tsv': self.create_file_path('proximal_variants_tsv'),
                'flanking_nucleotide_bases': self.flanking_bases * 4,
                'biotypes': self.biotypes,
                'allow_incomplete_transcripts': self.allow_incomplete_transcripts,
            }
            converter = VcfConverter(**params, pipeline_type='pVACseq')
            converter.execute()
            print("Completed")

    # creates transcripts.fa
    def tsv_to_fasta(self):
        if self.file_exists(self.create_file_path('fasta'), 'Variant fasta'):
            pass
        else:
            print('Creating variant fastas')
            params = {
                'input_file' : self.create_file_path('tsv'),
                'output_file': self.create_file_path('fasta'),
                'downstream_sequence_length': self.downstream_sequence_length,
                'proximal_variants_file': None if self.proximal_variants_vcf is None else self.create_file_path('proximal_variants_tsv')
            }
            variant_to_fasta = VariantToFasta(**params)
            variant_to_fasta.execute()
            print('Completed')

    def fasta_to_kmers(self):
        for el in self.choose_final_lengths():
            fasta_file = os.path.join(self.output_dir, f'{self.sample_name}.{el}.fa')
            if os.path.exists(fasta_file):
                print(f'{el}mer fasta already exists. Skipping.')
                continue
            else:
                print(f'Generating {el}mer peptides variant sequences')
                kmer_params = {
                    'fasta': self.create_file_path('fasta'),
                    'output_dir': self.output_dir,
                    'epitope_length': el,
                    'sample_name': self.sample_name,
                }
                fasta = VariantFastaToKmers(**kmer_params)
                fasta.execute()
                print('Completed')
