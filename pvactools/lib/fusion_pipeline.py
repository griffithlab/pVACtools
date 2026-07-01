import sys
import shutil
import os
import pandas as pd
from pvactools.lib.fusion_to_fasta import FusionToFasta
from pvactools.lib.fasta_to_kmers import FusionFastaToKmers
from pvactools.lib.combine_inputs import CombineInputs
from pvactools.lib.input_file_converter import FusionInputConverter

class FusionPipeline:
    def __init__(self, **kwargs):
        self.input_file = kwargs['input_file']
        self.sample_name = kwargs.pop('sample_name', "tmp")
        self.output_dir = kwargs['output_dir']
        self.transcript_fasta = kwargs['transcript_fasta']
        self.starfusion_file = kwargs.pop('starfusion_file', None)
        self.class_i_epitope_length = kwargs.pop('class_i_epitope_length', None)
        self.class_ii_epitope_length = kwargs.pop('class_ii_epitope_length', None)
        self.class_i_hla = kwargs.pop('class_i_hla', None)
        self.class_ii_hla = kwargs.pop('class_ii_hla', None)

    @staticmethod
    def file_exists(file_path: str, file_type: str):
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            print(f"{file_type} file already exists. Skipping.")
            exists = True
        else:
            exists = False
        return exists

    def execute(self):
        self.generate_fasta()
        self.fasta_to_kmers()

    def generate_fasta(self):
        self.input_to_tsv()
        self.fusion_to_fasta()

    def create_file_path(self, key):
        inputs = {
            'tsv': '.tsv',
            'fasta': '.transcripts.fa',
        }
        file_name = os.path.join(self.output_dir, self.sample_name + inputs[key])

        return file_name

    def input_to_tsv(self):
        if self.file_exists(self.create_file_path('tsv'), 'TSV'):
            pass
        else:
            print("Converting Fusion file to TSV")
            params = {
                'input_file' : self.input_file,
                'output_file': self.create_file_path('tsv'),
                'starfusion_file': self.starfusion_file
            }
            converter = FusionInputConverter(**params)
            converter.execute()
            print("Completed")

    # creates transcripts.fa
    def fusion_to_fasta(self):
        if self.file_exists(self.create_file_path('fasta'), 'Fusion fasta'):
            pass
        else:
            print('Creating fusion fastas')
            params = {
                'input_file' : self.create_file_path('tsv'),
                'transcript_fasta': self.transcript_fasta,
                'output_file': self.create_file_path('fasta'),
            }
            fusion_to_fasta = FusionToFasta(**params)
            fusion_to_fasta.execute()
            print('Completed')

    def fasta_to_kmers(self):
        for el in self.choose_final_lengths():
            fasta_file = os.path.join(self.output_dir, f'{self.sample_name}.{el}.fa')
            if os.path.exists(fasta_file):
                print(f'{el}mer fasta already exists. Skipping.')
                continue
            else:
                print(f'Generating {el}mer peptides fusion sequences')
                kmer_params = {
                    'fasta': self.create_file_path('fasta'),
                    'output_dir': self.output_dir,
                    'epitope_length': el,
                    'sample_name': self.sample_name,
                }
                fasta = FusionFastaToKmers(**kmer_params)
                fasta.execute()
                print('Completed')

    def choose_final_lengths(self):
        if not self.class_i_hla:
            lengths = self.class_ii_epitope_length
        elif not self.class_ii_hla:
            lengths = self.class_i_epitope_length
        else:
            lengths = self.class_i_epitope_length + self.class_ii_epitope_length
        return lengths
