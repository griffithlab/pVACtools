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

from pvactools.lib.fasta_generator import FusionFastaGenerator
from pvactools.lib.input_file_converter import FusionInputConverter
from pvactools.lib.calculate_manufacturability import CalculateManufacturability

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacfuse generate_protein_fasta",
        description="Generate an annotated fasta file from AGFusion or Arriba output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input",
        help="An AGFusion output directory or Arriba fusion.tsv output file."
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
        help = "A pVACfuse all_epitopes, filtered, or aggregated TSV file with epitopes to use for subsetting the input file to peptides of interest. Only the peptide sequences for the epitopes in the TSV will be used when creating the FASTA."
    )
    parser.add_argument(
        "--aggregate-report-evaluation",
        help="When running with an aggregate report input TSV, only include variants with this evaluation. Valid values for this field are Accept, Reject, Pending, and Review. Specifiy multiple values as a comma-separated list to include multiple evaluation states.",
        default='Accept',
        type=lambda s:[e for e in s.split(',')],
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshift fusion when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence."
    )
    return parser

def convert_fusion_input(input_file, temp_dir, starfusion_file):
    print("Converting Fusion file to TSV")
    tsv_file = os.path.join(temp_dir, 'tmp.tsv')
    convert_params = {
        'input_file' : input_file,
        'output_file': tsv_file,
        'starfusion_file': starfusion_file
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
        'downstream_sequence_length': downstream_sequence_length,
        'trim_invalid_characters'   : True,
    }
    fasta_generator = FusionFastaGenerator(**generate_fasta_params)
    fasta_generator.execute()
    print("Completed")
    if save_tsv_file:
        shutil.copy(tsv_file, "{}.tsv".format(args.output_file))

def parse_input_tsv(input_tsv):
    if input_tsv is None:
        return (None, None)
    indexes = []
    with open(input_tsv, 'r') as fh:
        reader = csv.DictReader(fh, delimiter = "\t")
        if 'Index' in reader.fieldnames:
            for line in reader:
                indexes.append(line[Index])
            file_type = 'full'
        else:
            for line in reader:
                indexes.append(line)
            file_type = 'aggregated'
    return (indexes, file_type)

def parse_files(output_file, temp_dir, input_tsv, aggregate_report_evaluation):
    print("Parsing the Variant Peptide FASTA and Key File")
    fasta_file_path = os.path.join(temp_dir, 'tmp.fasta')
    fasta_key_file_path = os.path.join(temp_dir, 'tmp.fasta.key')

    with open(fasta_key_file_path, 'r') as fasta_key_file:
        keys = yaml.load(fasta_key_file, Loader=yaml.FullLoader)

    (tsv_indexes, file_type) = parse_input_tsv(input_tsv)

    dataframe = OrderedDict()
    output_records = []
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        ids = keys[int(record.id)]
        for record_id in ids:
            if tsv_indexes is not None:
                if file_type == 'full':
                    if record_id not in tsv_indexes:
                        continue
                else:
                    matches = [r for r in tsv_indexes if r['ID'] == record_id and r['Evaluation'] in aggregate_report_evaluation]
                    if len(matches) == 0:
                        continue
            new_record = SeqRecord(record.seq, id=record_id, description=record_id)
            output_records.append(new_record)

    if tsv_indexes is not None:
        ordered_output_records = []
        for tsv_index in tsv_indexes:
            if file_type == 'full':
                records = [r for r in output_records if r.id == tsv_index]
            else:
                records = [r for r in output_records if r.id == tsv_index['ID']]
            ordered_output_records.extend(records)
        output_records = ordered_output_records

    SeqIO.write(output_records, output_file, "fasta")
    print("Completed")

def main(args_input = sys.argv[1:], save_tsv_file=False, starfusion_file=None):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.downstream_sequence_length == 'full':
        downstream_sequence_length = None
    elif args.downstream_sequence_length.isdigit():
        downstream_sequence_length = int(args.downstream_sequence_length)
    else:
        sys.exit("The downstream sequence length needs to be a positive integer or 'full'")

    temp_dir = tempfile.mkdtemp()
    convert_fusion_input(args.input, temp_dir, starfusion_file)
    generate_fasta(args, downstream_sequence_length, temp_dir, save_tsv_file)
    parse_files(args.output_file, temp_dir, args.input_tsv, args.aggregate_report_evaluation)
    shutil.rmtree(temp_dir, ignore_errors=True)
    manufacturability_file = "{}.manufacturability.tsv".format(args.output_file)
    print("Calculating Manufacturability Metrics")
    CalculateManufacturability(args.output_file, manufacturability_file, 'fasta').execute()
    print("Completed")

if __name__ == '__main__':
    main()
