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

from pvactools.lib.fusion_pipeline import FusionPipeline
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
        "ref_fasta",
        help="A reference CDS FASTA file. Note: this input should match the build and Ensembl version used to create the fusion annotations."
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
        help = "A pVACfuse all_epitopes, filtered, or aggregated TSV file with epitopes to use for subsetting the input file to peptides of interest. Only the peptide sequences for the variants in the TSV will be used when creating the FASTA. When running with an aggregated TSV, the sequences will be further narrowed down to only include variants with the selected --aggregate-report-evaluation."
    )
    parser.add_argument(
        "--mutant-only",
        help="Only output mutant peptide sequences",
        default=False,
        action='store_true',
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

def generate_fasta(args, temp_dir):
    params = {
        'input_file': args.input,
        'output_dir': temp_dir,
        'transcript_fasta': args.ref_fasta
    }
    pipeline = FusionPipeline(**params)
    pipeline.generate_fasta()

def parse_input_tsv(input_tsv):
    if input_tsv is None:
        return (None, None)
    indexes = []
    with open(input_tsv, 'r') as fh:
        reader = csv.DictReader(fh, delimiter = "\t")
        if 'ID' in reader.fieldnames:
            for line in reader:
                indexes.append(line)
            file_type = 'aggregated'
        else:
            for line in reader:
                indexes.append(line['Index'])
            file_type = 'full'
    return (indexes, file_type)

def trim_sequences(args, temp_dir, downstream_sequence_length):
    print("Trimming Variant Peptide FASTA")
    fasta_file_path = os.path.join(temp_dir, "tmp.transcripts.fa")
    trimmed_fasta_file_path = os.path.join(temp_dir, "tmp.transcripts.trimmed.fa")

    records = {}
    keys = set()
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        records[record.id] = str(record.seq)
        keys.add(record.id.split('.', 1)[1])

    output_records = []
    for key in sorted(keys, key=lambda x: int(x.split('.', 1)[0])):
        mt_seq = records[f"MT.{key}"]
        wt5_seq = records[f"WT5.{key}"]
        if f"WT3.{key}" in records:
            wt3_seq = records[f"WT3.{key}"]
        else:
            wt3_seq = None
        _, variant_type, position = key.rsplit('.', 2)
        position = int(position)
        start_position = position - args.flanking_sequence_length
        if start_position < 0:
            start_position = 0
        if variant_type == 'frameshift_fusion':
            if downstream_sequence_length is None:
                trimmed_mt_seq = mt_seq[start_position:]
            else:
                trimmed_mt_seq = mt_seq[start_position:(position + downstream_sequence_length)]
            output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
            if not args.mutant_only:
                trimmed_wt5_seq = wt5_seq[start_position:(position + args.flanking_sequence_length)]
                output_records.append(SeqRecord(Seq(trimmed_wt5_seq), id=f"WT5.{key}", description=""))
        else:
            end_position = position + args.flanking_sequence_length
            trimmed_mt_seq = mt_seq[start_position:end_position]
            output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
            if not args.mutant_only:
                trimmed_wt5_seq = wt5_seq[start_position:end_position]
                output_records.append(SeqRecord(Seq(trimmed_wt5_seq), id=f"WT5.{key}", description=""))
                wt3_position = len(wt3_seq) - len(mt_seq[position:])
                wt3_start_position = wt3_position - args.flanking_sequence_length
                if wt3_start_position < 0:
                    wt3_start_position = 0
                trimmed_wt3_seq = wt3_seq[wt3_start_position:(wt3_position + args.flanking_sequence_length)]
                output_records.append(SeqRecord(Seq(trimmed_wt3_seq), id=f"WT3.{key}", description=""))

    SeqIO.write(output_records, trimmed_fasta_file_path, "fasta")
    print("Completed")

def filter_fasta(args, temp_dir):
    trimmed_fasta_file_path = os.path.join(temp_dir, "tmp.transcripts.trimmed.fa")
    filtered_fasta_file_path = os.path.join(temp_dir, "tmp.transcripts.filtered.fa")

    if args.input_tsv is None:
        shutil.copy(trimmed_fasta_file_path, filtered_fasta_file_path)
    else:
        print("Filtering Variant Peptide FASTA")
        (tsv_indexes, file_type) = parse_input_tsv(args.input_tsv)

        output_records = []
        for record in SeqIO.parse(trimmed_fasta_file_path, "fasta"):
            record_id = record.id.split('.', 1)[1]
            if file_type == 'full':
                if record_id not in tsv_indexes:
                    continue
            else:
                matches = [r for r in tsv_indexes if r['ID'] == record_id and r['Evaluation'] in args.aggregate_report_evaluation]
                if len(matches) == 0:
                    continue
            new_record = SeqRecord(record.seq, id=record.id, description="")
            output_records.append(new_record)

        ordered_output_records = []
        for tsv_index in tsv_indexes:
            if file_type == 'full':
                records = [r for r in output_records if r.id.split('.', 1)[1] == tsv_index]
            else:
                records = [r for r in output_records if r.id.split('.', 1)[1] == tsv_index['ID']]
            ordered_output_records.extend(records)
        output_records = ordered_output_records

        SeqIO.write(output_records, filtered_fasta_file_path, "fasta")
        print("Completed")

    return(filtered_fasta_file_path)

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
    generate_fasta(args, temp_dir)
    trimmed_fasta = trim_sequences(args, temp_dir, downstream_sequence_length)
    filtered_fasta = filter_fasta(args, temp_dir)
    shutil.copy(filtered_fasta, args.output_file)
    shutil.rmtree(temp_dir, ignore_errors=True)
    manufacturability_file = "{}.manufacturability.tsv".format(args.output_file)
    print("Calculating Manufacturability Metrics")
    CalculateManufacturability(args.output_file, manufacturability_file, 'fasta').execute()
    print("Completed")

if __name__ == '__main__':
    main()
