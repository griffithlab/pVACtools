import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)
import argparse
import tempfile
import os
import yaml
import csv
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from lib.fasta_generator import *
from lib.input_file_converter import *
from lib.calculate_manufacturability import *

def define_parser():
    parser = argparse.ArgumentParser("pvacseq generate_protein_fasta", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "input_vcf",
        help="A VEP-annotated single-sample VCF containing transcript, Wildtype protein sequence, and Downstream protein sequence information."
    )
    parser.add_argument(
        "sample_name",
        help="Tumor sample name."
    )
    parser.add_argument(
        "peptide_sequence_length", type=int,
        help="Length of the peptide sequence to use when creating the FASTA.",
    )
    parser.add_argument(
        "output_file",
        help="The output fasta file."
    )
    parser.add_argument(
        "--input-tsv",
        help = "A pVACseq all_epitopes or filtered TSV file with epitopes to use for subsetting the input VCF to peptides of interest. Only the peptide sequences for the epitopes in the TSV will be used when creating the FASTA."
    )
    parser.add_argument(
        "--mutant-only",
        help="Only output mutant peptide sequences",
        default=False,
        action='store_true',
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence."
    )
    parser.add_argument(
        "-p", "--phased_proximal_variants_vcf",
        help="A VCF with phased proximal variant information. Must be gzipped and tabix indexed. (default: None)"
    )

    return parser

def convert_vcf(input_vcf, temp_dir, sample_name, phased_proximal_variants_vcf, peptide_sequence_length):
    print("Converting VCF to TSV")
    tsv_file = os.path.join(temp_dir, 'tmp.tsv')
    convert_params = {
        'input_file' : input_vcf,
        'output_file': tsv_file,
        'sample_name': sample_name
    }

    if phased_proximal_variants_vcf is not None:
        convert_params['proximal_variants_vcf'] = phased_proximal_variants_vcf
        proximal_variants_tsv = os.path.join(temp_dir, 'phased_proximal_variants_tmp.tsv')
        convert_params['proximal_variants_tsv'] = proximal_variants_tsv
        proximal_variants_file = proximal_variants_tsv
        convert_params['peptide_length'] = peptide_sequence_length # peptide_sequence_length

    converter = VcfConverter(**convert_params)
    converter.execute()
    print("Completed")


def generate_fasta(peptide_sequence_length, downstream_sequence_length, temp_dir, sample_name, phased_proximal_variants_vcf):
    print("Generating Variant Peptide FASTA and Key File")
    tsv_file = os.path.join(temp_dir, 'tmp.tsv')
    proximal_variants_tsv = None
    if phased_proximal_variants_vcf is not None:
        proximal_variants_tsv =  os.path.join(temp_dir, 'phased_proximal_variants_tmp.tsv')
    fasta_file = os.path.join(temp_dir, 'tmp.fasta')
    fasta_key_file = os.path.join(temp_dir, 'tmp.fasta.key')
    generate_fasta_params = {
        'input_file'                : tsv_file,
        'sample_name'               : sample_name,
        'peptide_sequence_length'   : peptide_sequence_length,
        'epitope_length'            : 0,
        'output_file'               : fasta_file,
        'output_key_file'           : fasta_key_file,
        'downstream_sequence_length': downstream_sequence_length,
        'proximal_variants_file'    : proximal_variants_tsv
    }
    fasta_generator = FastaGenerator(**generate_fasta_params)
    fasta_generator.execute()
    print("Completed")

def parse_input_tsv(input_tsv):
    if input_tsv is None:
        return None
    indexes = []
    with open(input_tsv, 'r') as fh:
        reader = csv.DictReader(fh, delimiter = "\t")
        for line in reader:
            consequence = line['Variant Type']
            if consequence == 'FS':
                amino_acid_change_position = "{}{}/{}".format(line['Protein Position'], line['Reference'], line['Variant'])
            else:
                amino_acid_change_position = "{}{}".format(line['Protein Position'], line['Mutation'])
            index = '%s.%s.%s.%s' % (line['Gene Name'], line['Transcript'], consequence, amino_acid_change_position)
            indexes.append(index)
    return indexes

def parse_files(output_file, temp_dir, mutant_only, input_tsv):
    print("Parsing the Variant Peptide FASTA and Key File")
    fasta_file_path = os.path.join(temp_dir, 'tmp.fasta')
    fasta_key_file_path = os.path.join(temp_dir, 'tmp.fasta.key')

    with open(fasta_key_file_path, 'r') as fasta_key_file:
        keys = yaml.load(fasta_key_file, Loader=yaml.FullLoader)

    tsv_indexes = parse_input_tsv(input_tsv)

    dataframe = OrderedDict()
    output_records = []
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        ids = keys[int(record.id)]
        for record_id in ids:
            if mutant_only and record_id.startswith('WT.'):
                continue
            if tsv_indexes is not None:
                sequence_type, count, index = record_id.split('.', 2)
                if index not in tsv_indexes:
                    continue
            new_record = SeqRecord(record.seq, id=record_id, description=record_id)
            output_records.append(new_record)

    SeqIO.write(output_records, output_file, "fasta")
    print("Completed")

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if "." in args.sample_name:
            sys.exit("Sample name cannot contain '.'")

    if args.downstream_sequence_length == 'full':
        downstream_sequence_length = None
    elif args.downstream_sequence_length.isdigit():
        downstream_sequence_length = int(args.downstream_sequence_length)
    else:
        sys.exit("The downstream sequence length needs to be a positive integer or 'full'")


    print(args.phased_proximal_variants_vcf)
    temp_dir = tempfile.mkdtemp()
    convert_vcf(args.input_vcf, temp_dir, args.sample_name, args.phased_proximal_variants_vcf, args.peptide_sequence_length)
    generate_fasta(args.peptide_sequence_length, downstream_sequence_length, temp_dir, args.sample_name, args.phased_proximal_variants_vcf)
    parse_files(args.output_file, temp_dir, args.mutant_only, args.input_tsv)
    manufacturability_file = "{}.manufacturability.tsv".format(args.output_file)
    print("Calculating Manufacturability Metrics")
    CalculateManufacturability(args.output_file, manufacturability_file, 'fasta').execute()
    print("Completed")

if __name__ == '__main__':
    main()
