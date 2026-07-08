import sys
import argparse
import tempfile
import os
import shutil
import yaml
import csv
import re
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pvactools.lib.junction_to_kmer_pipeline import JunctionToKmerPipeline
from pvactools.lib.calculate_manufacturability import CalculateManufacturability
from pvactools.lib.run_utils import *

def define_parser():
    parser = argparse.ArgumentParser(
        "pvacsplice generate_protein_fasta",
        description="Generate an annotated fasta file from a RegTools junctions output TSV file with protein sequences of mutations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_file",
        help="RegTools junctions output TSV file"
    )
    parser.add_argument(
        "flanking_sequence_length", type=int,
        help="Number of amino acids to add on each side of the splice site when creating the FASTA.",
    )
    parser.add_argument(
        "output_file",
        help="The output fasta file."
    )
    parser.add_argument(
        "annotated_vcf",
        help="A VEP-annotated single- or multi-sample VCF containing genotype and transcript information."
        + "The VCF ma be gzipped (requires tabix index)."
    )
    parser.add_argument(
        "ref_fasta",
        help="A reference FASTA file. Note: this input should be the same as the RegTools vcf input."
    )
    parser.add_argument(
        "gtf_file",
        help="A reference GTF file. Note: this input should be the same as the RegTools gtf input."
    )
    parser.add_argument(
        "--input-tsv",
        help = "A pVACsplice all_epitopes, filtered, or aggregated TSV file with epitopes to use for subsetting the inputs to peptides of interest. Only the peptide sequences for the epitopes in the TSV will be used when creating the FASTA. When running with an aggregated TSV, the sequences will be further narrowed down to only include variants with the selected --aggregate-report-evaluation."
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
        "-j", "--junction-score", type=int,
        help="Junction Coverage Cutoff. Only sites above this read depth cutoff will be considered.",
        default=10
    )
    parser.add_argument(
        "-v", "--variant-distance", type=int,
        help="Regulatory variants can lie inside or outside of splicing junction."
        + "Maximum distance window (upstream and downstream) for a variant outside the junction.",
        default=100
    )
    parser.add_argument(
        "--anchor-types", type=pvacsplice_anchors(),
        help="The anchor types of junctions to use. Multiple anchors can be specified using a comma-separated list."
        + "Choices: A, D, NDA, DA, N",
        default=['A', 'D', 'NDA'],
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
        "-s", "--sample-name",
        help="The name of the sample being processed. Required when processing a multi-sample VCF and must be a sample ID in the input VCF #CHROM header line."
    )
    return parser

def parse_input_tsv(input_tsv):
    if input_tsv is None:
        return (None, None)
    with open(input_tsv, 'r') as fh:
        reader = csv.DictReader(fh, delimiter = "\t")
        if 'Index' in reader.fieldnames:
            indexes = parse_full_input_tsv(reader)
            file_type = 'full'
        else:
            indexes = parse_aggregated_input_tsv(reader)
            file_type = 'aggregated'
    return (indexes, file_type)

def parse_full_input_tsv(reader):
    indexes = []
    for line in reader:
        indexes.append(line['Index'])
    return indexes

def parse_aggregated_input_tsv(reader):
    indexes = []
    for line in reader:
        indexes.append(line)
    return indexes

def generate_fasta(args, temp_dir):
    junction_arguments = {
        'input_file_type'                  : 'junctions',
        'junctions_dir'                    : temp_dir,
        'input_file'                       : args.input_file,
        'gtf_file'                         : args.gtf_file,
        'save_gtf'                         : False,
        'sample_name'                      : args.sample_name,
        'ref_fasta'                        : args.ref_fasta,
        'annotated_vcf'                    : args.annotated_vcf,
        'pass_only'                        : args.pass_only,
        'biotypes'                         : args.biotypes,
        'allow_incomplete_transcripts'     : args.allow_incomplete_transcripts,
        'junction_score'                   : args.junction_score,
        'variant_distance'                 : args.variant_distance,
        'anchor_types'                     : args.anchor_types,
        'normal_sample_name'               : None,
        'keep_tmp_files'                   : False,
        'class_i_epitope_length'           : [],
        'class_ii_epitope_length'          : [],
        'class_i_hla'                      : [],
        'class_ii_hla'                     : [],
    }

    pipeline = JunctionToKmerPipeline(**junction_arguments)
    pipeline.generate_fasta()

    return pipeline.create_file_path('fasta')

def trim_sequences(args, temp_dir, transcript_fasta):
    trimmed_transcript_fasta = os.path.join(temp_dir, f"{args.sample_name}.transcripts.trimmed.fa")

    records = {}
    keys = set()
    for record in SeqIO.parse(transcript_fasta, "fasta"):
        records[record.id] = record.seq
        keys.add(record.id.split('.', 1)[1])

    output_records = []
    for key in sorted(keys, key=lambda x: int(x.split('.', 1)[0])):
        mt_sequence = records[f"ALT.{key}"]
        wt_sequence = records[f"WT.{key}"]
        if mt_sequence in wt_sequence:
            continue
        _, frameshift_status = key.rsplit('.', 1)
        if frameshift_status == 'inframe_splice_site':
            final_mt_sequence, final_wt_sequence = get_mutated_peptide_with_flanking_sequence(wt_sequence, mt_sequence, min(args.flanking_sequence_length, len(wt_sequence)-1, len(mt_sequence)-1))
        elif frameshift_status == 'frameshift_splice_site':
            final_mt_sequence, final_wt_sequence = get_mutated_frameshift_peptide_with_flanking_sequence(wt_sequence, mt_sequence, min(args.flanking_sequence_length, len(wt_sequence)-1, len(mt_sequence)-1))
        else:
            raise Exception("Unexpected frameshift status {} for record {}. Skipping".format(frameshift_status, identifier))
        if final_mt_sequence and final_wt_sequence:
            output_records.append(SeqRecord(final_mt_sequence, id=f"MT.{key}", description=""))
            if not args.mutant_only:
                output_records.append(SeqRecord(final_wt_sequence, id=f"ALT.{key}", description=""))

    SeqIO.write(output_records, trimmed_transcript_fasta, "fasta")
    return trimmed_transcript_fasta

def filter_fasta(args, temp_dir):
    trimmed_transcript_fasta = os.path.join(temp_dir, f"{args.sample_name}.transcripts.trimmed.fa")
    filtered_transcript_fasta = os.path.join(temp_dir, f"{args.sample_name}.transcripts.filtered.fa")

    if args.input_tsv is None:
        shutil.copy(trimmed_transcript_fasta, filtered_transcript_fasta)
    else:
        (tsv_indexes, tsv_file_type) = parse_input_tsv(args.input_tsv)

        output_records = []
        for record in SeqIO.parse(trimmed_transcript_fasta, "fasta"):
            record_id = record.id.split('.', 1)[1]
            if tsv_file_type == 'full':
                if record_id not in tsv_indexes:
                    continue
            else:
                matches = [i for i in tsv_indexes if i['ID'] == record_id and i['Evaluation'] in args.aggregate_report_evaluation]
                if len(matches) == 0:
                    continue
            new_record = SeqRecord(record.seq, id=record.id, description="")
            output_records.append(new_record)

        ordered_output_records = []
        for tsv_index in tsv_indexes:
            if tsv_file_type == 'full':
                records = [r for r in output_records if r.id.split('.', 1)[1] == tsv_index]
            else:
                records = [r for r in output_records if r.id.split('.', 1)[1] == tsv_index['ID']]
            ordered_output_records.extend(records)
        output_records = ordered_output_records

        SeqIO.write(output_records, filtered_transcript_fasta, "fasta")

    return filtered_transcript_fasta

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if not (set(args.aggregate_report_evaluation)).issubset(set(['Accept', 'Reject', 'Review', 'Pending'])):
        sys.exit("Aggregate report evaluation ({}) contains invalid values.".format(args.aggregate_report_evaluation))

    temp_dir = tempfile.mkdtemp()
    transcript_fasta = generate_fasta(args, temp_dir)
    trimmed_fasta = trim_sequences(args, temp_dir, transcript_fasta)
    filtered_fasta = filter_fasta(args, temp_dir)
    shutil.copy(filtered_fasta, args.output_file)
    shutil.rmtree(temp_dir, ignore_errors=True)
    manufacturability_file = "{}.manufacturability.tsv".format(args.output_file)
    print("Calculating Manufacturability Metrics")
    CalculateManufacturability(args.output_file, manufacturability_file, 'fasta').execute()
    print("Completed")

if __name__ == '__main__':
    main()
