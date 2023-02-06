import os
import shutil
import sys
from pyfaidx import Fasta, FastaVariant

def create_personal_fasta(fasta_path, alt_fasta_path, annotated_vcf, sample_name):

    size1 = os.path.getsize(fasta_path)

    #if not os.path.exists(alt_fasta_path):
    print('Building personalized fasta.')
    shutil.copy(fasta_path, alt_fasta_path)
    personal_fasta = FastaVariant(alt_fasta_path, annotated_vcf, sample=sample_name)
    size2 = os.path.getsize(alt_fasta_path)
    if os.path.exists(alt_fasta_path) and size1 == size2:
        print('Completed')
    
    elif os.path.exists(alt_fasta_path) and size1 > os.path.getsize(alt_fasta_path):
        print('Personalized fasta file is incomplete. Trying again.')
        shutil.copy(fasta_path, alt_fasta_path)
        personal_fasta = FastaVariant(alt_fasta_path, annotated_vcf, sample=sample_name)
        size2 = os.path.getsize(alt_fasta_path)
        try: 
            size1 == size2
            print('Completed')
        except:
            sys.exit('Error: Reference and personalized fasta files are different sizes.')

    # elif os.path.exists(alt_fasta_path) and size1 == os.path.getsize(alt_fasta_path):
    #     print('Personalized fasta already exists. Skipping.')
    #     personal_fasta = Fasta(alt_fasta_path)

    return personal_fasta

    
