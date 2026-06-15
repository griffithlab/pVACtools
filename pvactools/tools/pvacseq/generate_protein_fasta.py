import sys
import argparse
import tempfile
import os
import shutil
import yaml
import csv
import re
import json
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pvactools.lib.variant_pipeline import VariantPipeline
from pvactools.lib.calculate_manufacturability import CalculateManufacturability
from pvactools.lib.run_utils import *

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacseq generate_protein_fasta",
        description="Generate an annotated fasta file from a VCF with protein sequences of mutations and matching wildtypes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input_vcf",
        help="A VEP-annotated single- or multi-sample VCF containing genotype, transcript, "
            +"Wildtype protein sequence, and Frameshift protein sequence information."
            +"The VCF may be gzipped (requires tabix index)."
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
        help = "A pVACseq all_epitopes, filtered, or aggregated TSV file with epitopes to use for subsetting the input VCF to peptides of interest. Only the peptide sequences for the epitopes in the TSV will be used when creating the FASTA. When running with an aggregated TSV, the sequences will be further narrowed down to only include variants with the selected --aggregate-report-evaluation."
    )
    parser.add_argument(
        "-p", "--phased-proximal-variants-vcf",
        help="A VCF with phased proximal variant information to incorporate into the predicted fasta sequences. Must be gzipped and tabix indexed."
    )
    parser.add_argument(
        '--pass-only',
        help="Only process VCF entries with a PASS status.",
        default=False,
        action='store_true',
    )
    parser.add_argument(
        "--biotypes", type=lambda s:[a for a in s.split(',')],
        help="A list of biotypes to use for pre-filtering transcripts for processing in the pipeline.",
        default=['protein_coding']
    )
    parser.add_argument(
        "--allow-incomplete-transcripts",
        help="By default, transcripts annotated with incomplete CDS (i.e., 'cds_start_NF' or 'cds_end_NF' flags in the VEP CSQ field) "
                + "are excluded from analysis, as they often produce invalid protein sequences. "
                + "Use this flag to allow candidates from such transcripts. Only peptides that do not contain 'X' will be included. "
                + "These candidates will be deprioritized relative to those from transcripts without incomplete CDS flags.",
        default=False,
        action='store_true'
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
        type=aggregate_report_evaluations(),
    )
    parser.add_argument(
        "-d", "--downstream-sequence-length",
        default="1000",
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence."
    )
    parser.add_argument(
        "-s", "--sample-name",
        help="The name of the sample being processed. Required when processing a multi-sample VCF and must be a sample ID in the input VCF #CHROM header line."
    )
    return parser

def generate_fasta(args, temp_dir, downstream_sequence_length):
    params = {
        'output_dir'                  : temp_dir,
        'input_file'                  : args.input_vcf,
        'sample_name'                 : args.sample_name or 'tmp',
        'pass_only'                   : args.pass_only,
        'proximal_variants_vcf'       : args.phased_proximal_variants_vcf,
        'biotypes'                    : args.biotypes,
        'allow_incomplete_transcripts': args.allow_incomplete_transcripts,
        'downstream_sequence_length'  : downstream_sequence_length,
        'flanking_bases'              : args.flanking_sequence_length,
    }
    pipeline = VariantPipeline(**params)
    pipeline.generate_fasta()

def parse_input_tsv(input_tsv):
    if input_tsv is None:
        return (None, None)
    indexes = []
    with open(input_tsv, 'r') as fh:
        reader = csv.DictReader(fh, delimiter = "\t")
        if 'Best Peptide' in reader.fieldnames:
            for line in reader:
                indexes.append(line)
            file_type = 'aggregated'
        else:
            for line in reader:
                indexes.append(line['Index'])
            file_type = 'full'
    return (indexes, file_type)

def trim_sequences(args=None, temp_dir=None, fasta_file_path=None, trimmed_fasta_file_path=None, flanking_sequence_length=None, mutant_only=None):
    print("Trimming Variant Peptide FASTA")
    if fasta_file_path is None:
        fasta_file_path = os.path.join(temp_dir, f"{args.sample_name or 'tmp'}.transcripts.fa")
    if trimmed_fasta_file_path is None:
        trimmed_fasta_file_path = os.path.join(temp_dir, f"{args.sample_name or 'tmp'}.transcripts.trimmed.fa")
    if flanking_sequence_length is None:
        flanking_sequence_length = args.flanking_sequence_length
    if mutant_only is None:
        mutant_only = args.mutant_only

    records = {}
    keys = set()
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        records[record.id] = str(record.seq)
        keys.add(record.id.split('.', 1)[1])

    output_records = []
    for key in sorted(keys, key=lambda x: int(x.split('.', 1)[0])):
        mt_seq = records[f"MT.{key}"]
        wt_seq = records[f"WT.{key}"]
        _, variant_type, aa_change = key.rsplit('.', 2)
        position = int(re.split('[A-Z|-]', aa_change)[0])
        start_position = position - flanking_sequence_length - 1
        if start_position < 0:
            start_position = 0
        end_position = position + flanking_sequence_length
        if variant_type == 'missense':
            trimmed_mt_seq = mt_seq[start_position:end_position]
            output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
            if not mutant_only:
                trimmed_wt_seq = wt_seq[start_position:end_position]
                output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))
        elif variant_type == 'FS':
            trimmed_mt_seq = mt_seq[start_position:]
            output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
            if not mutant_only:
                trimmed_wt_seq = wt_seq[start_position:end_position]
                output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))
        elif variant_type == 'inframe_del':
            trimmed_mt_seq = mt_seq[start_position:end_position]
            output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
            if not mutant_only:
                offset = len(wt_seq) - len(mt_seq)
                trimmed_wt_seq = wt_seq[start_position:(end_position + offset)]
                output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))
        else:
            offset = len(mt_seq) - len(wt_seq)
            trimmed_mt_seq = mt_seq[start_position:(end_position + offset)]
            output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
            if not mutant_only:
                trimmed_wt_seq = wt_seq[start_position:end_position]
                output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))

    SeqIO.write(output_records, trimmed_fasta_file_path, "fasta")
    print("Completed")

def filter_fasta(args, temp_dir):
    trimmed_fasta_file_path = os.path.join(temp_dir, f"{args.sample_name or 'tmp'}.transcripts.trimmed.fa")
    filtered_fasta_file_path = os.path.join(temp_dir, f"{args.sample_name or 'tmp'}.transcripts.filtered.fa")

    if args.input_tsv is None:
        shutil.copy(trimmed_fasta_file_path, filtered_fasta_file_path)
    else:
        print("Filtering Variant Peptide FASTA")
        (tsv_indexes, file_type) = parse_input_tsv(args.input_tsv)

        output_records = []
        for record in SeqIO.parse(trimmed_fasta_file_path, "fasta"):
            record_id = record.id.split('.', 1)[1]
            description = ""
            if file_type == 'full':
                if record_id not in tsv_indexes:
                    continue
            else:
                matches = [r for r in tsv_indexes if r['Index'] == record_id and r['Evaluation'] in args.aggregate_report_evaluation]
                if len(matches) == 0:
                    continue
                elif len(matches) == 1 and record.id.startswith('MT.'):
                    description = json.dumps({ 'Best Peptide': matches[0]['Best Peptide'] })
                elif len(matches) > 1:
                    import pdb
                    pdb.set_trace()
            new_record = SeqRecord(record.seq, id=record.id, description=description)
            output_records.append(new_record)

        ordered_output_records = []
        for tsv_index in tsv_indexes:
            if file_type == 'full':
                records = [r for r in output_records if r.id.split('.', 1)[1] == tsv_index]
            else:
                records = [r for r in output_records if r.id.split('.', 1)[1] == tsv_index['Index']]
            ordered_output_records.extend(records)
        output_records = ordered_output_records

        SeqIO.write(output_records, filtered_fasta_file_path, "fasta")
        print("Completed")

    return(filtered_fasta_file_path)

#todo: remove after figuring out how to incorporate order sheet creation
def run_generate_protein_fasta(
    input_vcf,
    flanking_sequence_length,
    output_file,
    input_tsv=None,
    phased_proximal_variants_vcf=None,
    pass_only=False,
    biotypes=['protein_coding'],
    allow_incomplete_transcripts=False,
    mutant_only=False,
    aggregate_report_evaluation=['Accept'],
    downstream_sequence_length="1000",
    sample_name=None,
    peptide_ordering_form=False
):

    return
    try:
        generate_fasta(
            downstream_sequence_length=downstream_sequence_length,
            temp_dir=temp_dir,
        )

        parse_files(
            output_file=output_file,
            temp_dir=temp_dir,
            mutant_only=mutant_only,
            input_tsv=input_tsv,
            aggregate_report_evaluation=aggregate_report_evaluation
        )

        if peptide_ordering_form:
            parse_files(
                output_file=f"{output_file}_combined",
                temp_dir=temp_dir,
                mutant_only=not mutant_only,
                input_tsv=input_tsv,
                aggregate_report_evaluation=aggregate_report_evaluation
            )

        manufacturability_file = f"{output_file}.manufacturability.tsv"
        print("Calculating Manufacturability Metrics")
        CalculateManufacturability(output_file, manufacturability_file, 'fasta').execute()
        print("Completed")

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

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
    generate_fasta(
        args,
        temp_dir,
        downstream_sequence_length,
    )
    trimmed_fasta = trim_sequences(args, temp_dir)
    filtered_fasta = filter_fasta(args, temp_dir)
    shutil.copy(filtered_fasta, args.output_file)
    shutil.rmtree(temp_dir, ignore_errors=True)
    manufacturability_file = "{}.manufacturability.tsv".format(args.output_file)
    print("Calculating Manufacturability Metrics")
    CalculateManufacturability(args.output_file, manufacturability_file, 'fasta').execute()
    print("Completed")

if __name__ == '__main__':
    main()
