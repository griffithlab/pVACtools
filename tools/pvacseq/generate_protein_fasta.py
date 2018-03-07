import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)
import argparse
import tempfile
import os
import yaml
from collections import OrderedDict
from lib.fasta_generator import *
from lib.input_file_converter import *

def define_parser():
    parser = argparse.ArgumentParser("pvacseq generate_protein_fasta")

    parser.add_argument(
        "input_file",
        help="A VEP-annotated single-sample VCF containing transcript, Wildtype protein sequence, and Downstream protein sequence information"
    )
    parser.add_argument(
        "peptide_sequence_length", type=int,
        help="Length of the peptide sequence to use when creating the FASTA.",
    )
    parser.add_argument(
        "output_file",
        help="The output fasta file"
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence. Default: 1000"
    )
    return parser

def convert_vcf(input_file, temp_dir):
    print("Converting VCF to TSV")
    tsv_file = os.path.join(temp_dir, 'tmp.tsv')
    convert_params = {
        'input_file' : input_file,
        'output_file': tsv_file,
        'gene_expn_file'             : None,
        'transcript_expn_file'       : None,
        'normal_snvs_coverage_file'  : None,
        'normal_indels_coverage_file': None,
        'tdna_snvs_coverage_file'    : None,
        'tdna_indels_coverage_file'  : None,
        'trna_snvs_coverage_file'    : None,
        'trna_indels_coverage_file'  : None,
    }
    converter = VcfConverter(**convert_params)
    converter.execute()
    print("Completed")

def generate_fasta(peptide_sequence_length, downstream_sequence_length, temp_dir):
    print("Generating Variant Peptide FASTA and Key File")
    tsv_file = os.path.join(temp_dir, 'tmp.tsv')
    fasta_file = os.path.join(temp_dir, 'tmp.fasta')
    fasta_key_file = os.path.join(temp_dir, 'tmp.fasta.key')
    generate_fasta_params = {
        'input_file'                : tsv_file,
        'peptide_sequence_length'   : peptide_sequence_length,
        'epitope_length'            : 0,
        'output_file'               : fasta_file,
        'output_key_file'           : fasta_key_file,
        'downstream_sequence_length': downstream_sequence_length
    }
    fasta_generator = FastaGenerator(**generate_fasta_params)
    fasta_generator.execute()
    print("Completed")

def parse_files(output_file, temp_dir):
    print("Parsing the Variant Peptide FASTA and Key File")
    fasta_file_path = os.path.join(temp_dir, 'tmp.fasta')
    fasta_key_file_path = os.path.join(temp_dir, 'tmp.fasta.key')

    with open(fasta_key_file_path, 'r') as fasta_key_file:
        keys = yaml.load(fasta_key_file)

    dataframe = OrderedDict()
    with open(fasta_file_path, 'r') as fasta_file:
        for line in fasta_file:
            key      = line.rstrip().replace(">","")
            sequence = fasta_file.readline().rstrip()
            ids      = keys[int(key)]
            for id in ids:
                (type, index) = id.split('.', 1)
                if index not in dataframe:
                    dataframe[index] = {}
                dataframe[index][type] = sequence

    with open(output_file, 'w') as parsed_fasta_file:
        for index, sequences in dataframe.items():
            for type in ('WT', 'MT'):
                parsed_fasta_file.write(">%s.%s\n" % (type, index))
                parsed_fasta_file.write("%s\n" % sequences[type])
    print("Completed")

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.downstream_sequence_length == 'full':
        downstream_sequence_length = None
    elif args.downstream_sequence_length.isdigit():
        downstream_sequence_length = int(args.downstream_sequence_length)
    else:
        sys.exit("The downstream sequence length needs to be a positive integer or 'full'")

    temp_dir = tempfile.mkdtemp()
    convert_vcf(args.input_file, temp_dir)
    generate_fasta(args.peptide_sequence_length, downstream_sequence_length, temp_dir)
    parse_files(args.output_file, temp_dir)

if __name__ == '__main__':
    main()
