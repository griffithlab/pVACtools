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
from Bio.Alphabet import IUPAC

from pvactools.lib.fasta_generator import FastaGenerator
from pvactools.lib.input_file_converter import VcfConverter
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
        help = "A pVACseq all_epitopes, filtered, or aggregated TSV file with epitopes to use for subsetting the input VCF to peptides of interest. Only the peptide sequences for the epitopes in the TSV will be used when creating the FASTA."
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
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file. "
            + "Use 'full' to include the full downstream sequence."
    )
    parser.add_argument(
        "-s", "--sample-name",
        help="The name of the sample being processed. Required when processing a multi-sample VCF and must be a sample ID in the input VCF #CHROM header line."
    )
    return parser

def convert_vcf(input_vcf, temp_dir, sample_name, phased_proximal_variants_vcf, flanking_sequence_length, pass_only):
    print("Converting VCF to TSV")
    tsv_file = os.path.join(temp_dir, 'tmp.tsv')
    convert_params = {
        'input_file' : input_vcf,
        'output_file': tsv_file,
    }
    if sample_name is not None:
        convert_params['sample_name'] = sample_name
    if phased_proximal_variants_vcf is not None:
        convert_params['proximal_variants_vcf'] = phased_proximal_variants_vcf
        proximal_variants_tsv = os.path.join(temp_dir, 'proximal_variants.tsv')
        convert_params['proximal_variants_tsv'] = proximal_variants_tsv
        convert_params['flanking_bases'] = (flanking_sequence_length + 1)* 4
    else:
        proximal_variants_tsv = None
    if pass_only:
        convert_params['pass_only'] = pass_only
    converter = VcfConverter(**convert_params)
    converter.execute()
    print("Completed")
    return proximal_variants_tsv

def generate_fasta(flanking_sequence_length, downstream_sequence_length, temp_dir, proximal_variants_tsv):
    print("Generating Variant Peptide FASTA and Key File")
    tsv_file = os.path.join(temp_dir, 'tmp.tsv')
    fasta_file = os.path.join(temp_dir, 'tmp.fasta')
    fasta_key_file = os.path.join(temp_dir, 'tmp.fasta.key')
    generate_fasta_params = {
        'input_file'                : tsv_file,
        'flanking_sequence_length'  : flanking_sequence_length,
        'epitope_length'            : 0,
        'output_file'               : fasta_file,
        'output_key_file'           : fasta_key_file,
        'downstream_sequence_length': downstream_sequence_length,
        'proximal_variants_file'    : proximal_variants_tsv,
    }
    fasta_generator = FastaGenerator(**generate_fasta_params)
    fasta_generator.execute()
    print("Completed")

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

def parse_files(output_file, temp_dir, mutant_only, input_tsv, aggregate_report_evaluation):
    print("Parsing the Variant Peptide FASTA and Key File")
    fasta_file_path = os.path.join(temp_dir, 'tmp.fasta')
    fasta_key_file_path = os.path.join(temp_dir, 'tmp.fasta.key')

    with open(fasta_key_file_path, 'r') as fasta_key_file:
        keys = yaml.load(fasta_key_file, Loader=yaml.FullLoader)

    (tsv_indexes, tsv_file_type) = parse_input_tsv(input_tsv)

    dataframe = OrderedDict()
    output_records = []
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        ids = keys[int(record.id)]
        for record_id in ids:
            if mutant_only and record_id.startswith('WT.'):
                continue
            if tsv_indexes is not None:
                if tsv_file_type == 'full':
                    sequence_type, index = record_id.split('.', 1)
                    if index not in tsv_indexes:
                        continue
                else:
                    (rest_record_id, variant_type, aa_change) = record_id.rsplit(".", 2)
                    (peptide_type, count, gene, transcript) = rest_record_id.split(".", 3)
                    (parsed_aa_change, _, _, _) = index_to_aggregate_report_aa_change(aa_change, variant_type)
                    matches = [i for i in tsv_indexes if i['Best Transcript'] == transcript and i['AA Change'] == parsed_aa_change and i['Evaluation'] in aggregate_report_evaluation]
                    if len(matches) == 0:
                        continue
            new_record = SeqRecord(record.seq, id=record_id, description=record_id)
            output_records.append(new_record)

    SeqIO.write(output_records, output_file, "fasta")
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

    if not (set(args.aggregate_report_evaluation)).issubset(set(['Accept', 'Reject', 'Review', 'Pending'])):
        sys.exit("Aggregate report evaluation ({}) contains invalid values.".format(args.aggregate_report_evaluation))

    temp_dir = tempfile.mkdtemp()
    proximal_variants_tsv = convert_vcf(args.input_vcf, temp_dir, args.sample_name, args.phased_proximal_variants_vcf, args.flanking_sequence_length, args.pass_only)
    generate_fasta(args.flanking_sequence_length, downstream_sequence_length, temp_dir, proximal_variants_tsv)
    parse_files(args.output_file, temp_dir, args.mutant_only, args.input_tsv, args.aggregate_report_evaluation)
    shutil.rmtree(temp_dir, ignore_errors=True)
    manufacturability_file = "{}.manufacturability.tsv".format(args.output_file)
    print("Calculating Manufacturability Metrics")
    CalculateManufacturability(args.output_file, manufacturability_file, 'fasta').execute()
    print("Completed")

if __name__ == '__main__':
    main()
