import sys
import argparse
import tempfile
import os
import shutil
import yaml
import csv
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from pvactools.lib.fasta_generator import FusionFastaGenerator
from pvactools.lib.input_file_converter import FusionInputConverter
from pvactools.lib.calculate_manufacturability import CalculateManufacturability

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacfuse generate_protein_fasta",
        description="Generate an annotated fasta file from AGFusion output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input",
        help="An AGfusion output directory."
    )
    parser.add_argument(
        "flanking_sequence_length", type=int,
        help="Number of amino acids to add on each side of the mutation when creating the FASTA.",
    )
    parser.add_argument(
        "output_file",
        help="The output fasta file."
    )
    parser.add_argument(
        "--input-tsv",
        help = "A pVACfuse all_epitopes or filtered TSV file with epitopes to use for subsetting the input file to peptides of interest. Only the peptide sequences for the epitopes in the TSV will be used when creating the FASTA."
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshift fusion when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence."
    )
    return parser

def convert_fusion_input(input_file, temp_dir):
    print("Converting Fusion file to TSV")
    tsv_file = os.path.join(temp_dir, 'tmp.tsv')
    convert_params = {
        'input_file' : input_file,
        'output_file': tsv_file,
    }
    converter = FusionInputConverter(**convert_params)
    converter.execute()
    print("Completed")

def generate_fasta(args, downstream_sequence_length, temp_dir, save_tsv_file):
    print("Generating Variant Peptide FASTA and Key File")
    tsv_file = os.path.join(temp_dir, 'tmp.tsv')
    fasta_file = os.path.join(temp_dir, 'tmp.fasta')
    fasta_key_file = os.path.join(temp_dir, 'tmp.fasta.key')
    generate_fasta_params = {
        'input_file'                : tsv_file,
        'flanking_sequence_length'  : args.flanking_sequence_length,
        'epitope_length'            : 0,
        'output_file'               : fasta_file,
        'output_key_file'           : fasta_key_file,
        'downstream_sequence_length': downstream_sequence_length
    }
    fasta_generator = FusionFastaGenerator(**generate_fasta_params)
    fasta_generator.execute()
    print("Completed")
    if save_tsv_file:
        shutil.copy(tsv_file, "{}.tsv".format(args.output_file))

def parse_input_tsv(input_tsv):
    if input_tsv is None:
        return None
    indexes = []
    with open(input_tsv, 'r') as fh:
        reader = csv.DictReader(fh, delimiter = "\t")
        for line in reader:
            indexes.append(line['Index'])
    return indexes

def parse_files(output_file, temp_dir, input_tsv):
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
            if tsv_indexes is not None:
                if record_id not in tsv_indexes:
                    continue
            new_record = SeqRecord(record.seq, id=record_id, description=record_id)
            output_records.append(new_record)

    SeqIO.write(output_records, output_file, "fasta")
    print("Completed")

def main(args_input = sys.argv[1:], save_tsv_file=False):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.downstream_sequence_length == 'full':
        downstream_sequence_length = None
    elif args.downstream_sequence_length.isdigit():
        downstream_sequence_length = int(args.downstream_sequence_length)
    else:
        sys.exit("The downstream sequence length needs to be a positive integer or 'full'")

    temp_dir = tempfile.mkdtemp()
    convert_fusion_input(args.input, temp_dir)
    generate_fasta(args, downstream_sequence_length, temp_dir, save_tsv_file)
    parse_files(args.output_file, temp_dir, args.input_tsv)
    shutil.rmtree(temp_dir, ignore_errors=True)
    manufacturability_file = "{}.manufacturability.tsv".format(args.output_file)
    print("Calculating Manufacturability Metrics")
    CalculateManufacturability(args.output_file, manufacturability_file, 'fasta').execute()
    print("Completed")

if __name__ == '__main__':
    main()
