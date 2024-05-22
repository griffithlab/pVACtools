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

from pvactools.lib.splice_pipeline import JunctionPipeline
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
        help = "A pVACsplice all_epitopes, filtered, or aggregated TSV file with epitopes to use for subsetting the inputs to peptides of interest. Only the peptide sequences for the epitopes in the TSV will be used when creating the FASTA."
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
        "--anchor-types", nargs="*",
        help="The anchor types of junctions to use. Multiple anchors can be specified using a comma-separated list."
        + "Choices: A, D, NDA, DA, N",
        default=['A', 'D', 'NDA'],
        choices=['A', 'D', 'NDA', 'DA', 'N']
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

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if not (set(args.aggregate_report_evaluation)).issubset(set(['Accept', 'Reject', 'Review', 'Pending'])):
        sys.exit("Aggregate report evaluation ({}) contains invalid values.".format(args.aggregate_report_evaluation))

    temp_dir = tempfile.mkdtemp()

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

    pipeline = JunctionPipeline(**junction_arguments)
    pipeline.vcf_to_tsv()
    pipeline.junction_to_fasta()

    transcript_fasta = pipeline.create_file_path('fasta')
    mt_sequences = {}
    wt_sequences = {}
    for record in SeqIO.parse(transcript_fasta, "fasta"):
        if record.id.startswith('WT.'):
            wt_sequences[record.id.split('.', 1)[1]] = record.seq
        if record.id.startswith('ALT.'):
            mt_sequences[record.id.split('.', 1)[1]] = record.seq

    final_sequences = {}
    for (index, mt_sequence) in mt_sequences.items():
        wt_sequence = wt_sequences[index]
        final_sequences[index] = get_mutated_peptide_with_flanking_sequence(wt_sequence, mt_sequence, args.flanking_sequence_length)

    (tsv_indexes, tsv_file_type) = parse_input_tsv(args.input_tsv)
    output_records = []
    for (index, sequence) in final_sequences.items():
        if tsv_indexes is not None:
            if tsv_file_type == 'full':
                if index not in tsv_indexes:
                    continue
            else:
                matches = [i for i in tsv_indexes if i['ID'] == index and i['Evaluation'] in args.aggregate_report_evaluation]
                if len(matches) == 0:
                    continue
        new_record = SeqRecord(sequence, id=index, description=index)
        output_records.append(new_record)

    SeqIO.write(output_records, args.output_file, "fasta")
    print("Completed")

    shutil.rmtree(temp_dir, ignore_errors=True)
    manufacturability_file = "{}.manufacturability.tsv".format(args.output_file)
    print("Calculating Manufacturability Metrics")
    CalculateManufacturability(args.output_file, manufacturability_file, 'fasta').execute()
    print("Completed")

if __name__ == '__main__':
    main()
