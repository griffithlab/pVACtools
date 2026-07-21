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

from pvactools.lib.variant_to_kmer_pipeline import VariantToKmerPipeline
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

class PvacseqGenerateProteinFasta():
    def __init__(self, **kwargs):
        self.input_vcf = kwargs.pop('input_vcf', None)
        self.sample_name = kwargs.pop('sample_name', 'tmp')
        if self.sample_name is None:
            self.sample_name = 'tmp'
        self.pass_only = kwargs.pop('pass_only', False)
        self.phased_proximal_variants_vcf = kwargs.pop('phased_proximal_variants_vcf', None)
        self.biotypes = kwargs.pop('biotypes', ['protein_coding'])
        self.allow_incomplete_transcripts = kwargs.pop('allow_incomplete_transcripts', False)
        self.downstream_sequence_length = kwargs.pop('downstream_sequence_length', 1000)
        self.flanking_sequence_length = kwargs.pop('flanking_sequence_length')
        self.mutant_only = kwargs.pop('mutant_only', False)
        self.aggregate_report_evaluation = kwargs.pop('aggregate_report_evaluation', ['Accept'])
        self.temp_dir = tempfile.mkdtemp()
        self.fasta_file_path = kwargs.pop('fasta_file_path', os.path.join(self.temp_dir, f"{self.sample_name}.transcripts.fa"))
        self.trimmed_fasta_file_path = kwargs.pop('trimmed_fasta_file_path', os.path.join(self.temp_dir, f"{self.sample_name}.transcripts.trimmed.fa"))
        self.filtered_fasta_file_path = os.path.join(self.temp_dir, f"{self.sample_name}.transcripts.filtered.fa")
        self.input_tsv = kwargs.pop('input_tsv', None)
        self.output_file = kwargs.pop('output_file', None)

    def generate_fasta(self):
        params = {
            'output_dir'                  : self.temp_dir,
            'input_file'                  : self.input_vcf,
            'sample_name'                 : self.sample_name,
            'pass_only'                   : self.pass_only,
            'proximal_variants_vcf'       : self.phased_proximal_variants_vcf,
            'biotypes'                    : self.biotypes,
            'allow_incomplete_transcripts': self.allow_incomplete_transcripts,
            'downstream_sequence_length'  : self.downstream_sequence_length,
            'flanking_bases'              : self.flanking_sequence_length,
        }
        pipeline = VariantToKmerPipeline(**params)
        pipeline.generate_fasta()

    def parse_input_tsv(self):
        if self.input_tsv is None:
            return (None, None)
        indexes = []
        with open(self.input_tsv, 'r') as fh:
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

    def trim_sequences(self):
        print("Trimming Variant Peptide FASTA")
        records = {}
        keys = set()
        for record in SeqIO.parse(self.fasta_file_path, "fasta"):
            records[record.id] = str(record.seq)
            keys.add(record.id.split('.', 1)[1])

        output_records = []
        for key in sorted(keys, key=lambda x: int(x.split('.', 1)[0])):
            mt_seq = records[f"MT.{key}"]
            wt_seq = records[f"WT.{key}"]
            _, variant_type, aa_change = key.rsplit('.', 2)
            position = int(re.split('[A-Z|-]', aa_change)[0])
            start_position = position - self.flanking_sequence_length - 1
            end_position = position + self.flanking_sequence_length
            if variant_type == 'missense':
                if start_position < 0:
                    start_position = 0
                trimmed_mt_seq = mt_seq[start_position:end_position]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    trimmed_wt_seq = wt_seq[start_position:end_position]
                    output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))
            elif variant_type == 'FS':
                if start_position < 0:
                    start_position = 0
                trimmed_mt_seq = mt_seq[start_position:]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    trimmed_wt_seq = wt_seq[start_position:end_position]
                    output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))
            elif variant_type == 'inframe_del':
                match = re.match(r"\d+(?:-\d+)?([A-Z]+)/([A-Z]+|-)", aa_change)
                wt_aa, mt_aa = match.groups()
                if wt_aa.startswith(mt_aa):
                    start_position = start_position + 1
                else:
                    end_position = end_position - 1
                if start_position < 0:
                    start_position = 0
                trimmed_mt_seq = mt_seq[start_position:(end_position)]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    offset = len(wt_seq) - len(mt_seq)
                    trimmed_wt_seq = wt_seq[start_position:(end_position + offset)]
                    output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))
            else:
                match = re.match(r"\d+(?:-\d+)?([A-Z]+|-)/([A-Z]+)", aa_change)
                wt_aa, mt_aa = match.groups()
                if mt_aa.startswith(wt_aa):
                    start_position = start_position + 1
                else:
                    end_position = end_position - 1
                if start_position < 0:
                    start_position = 0
                offset = len(mt_seq) - len(wt_seq)
                trimmed_mt_seq = mt_seq[start_position:(end_position + offset)]
                output_records.append(SeqRecord(Seq(trimmed_mt_seq), id=f"MT.{key}", description=""))
                if not self.mutant_only:
                    trimmed_wt_seq = wt_seq[start_position:(end_position)]
                    output_records.append(SeqRecord(Seq(trimmed_wt_seq), id=f"WT.{key}", description=""))

        SeqIO.write(output_records, self.trimmed_fasta_file_path, "fasta")
        print("Completed")

    def filter_fasta(self):
        if self.input_tsv is None:
            shutil.copy(self.trimmed_fasta_file_path, self.filtered_fasta_file_path)
        else:
            print("Filtering Variant Peptide FASTA")
            (tsv_indexes, file_type) = self.parse_input_tsv()

            output_records = []
            for record in SeqIO.parse(self.trimmed_fasta_file_path, "fasta"):
                record_id = record.id.split('.', 1)[1]
                description = ""
                if file_type == 'full':
                    if record_id not in tsv_indexes:
                        continue
                else:
                    matches = [r for r in tsv_indexes if r['Index'] == record_id and r['Evaluation'] in self.aggregate_report_evaluation]
                    if len(matches) == 0:
                        continue
                    elif len(matches) > 0 and record.id.startswith('MT.'):
                        description = json.dumps({ 'Best Peptide': matches[0]['Best Peptide'] })
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

            SeqIO.write(output_records, self.filtered_fasta_file_path, "fasta")
            print("Completed")

    def execute(self):
        self.generate_fasta()
        self.trim_sequences()
        self.filter_fasta()
        shutil.copy(self.filtered_fasta_file_path, self.output_file)
        shutil.rmtree(self.temp_dir, ignore_errors=True)
        manufacturability_file = "{}.manufacturability.tsv".format(self.output_file)
        print("Calculating Manufacturability Metrics")
        CalculateManufacturability(self.output_file, manufacturability_file, 'fasta').execute()
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

    params = {
        'input_vcf': args.input_vcf,
        'sample_name': args.sample_name,
        'pass_only': args.pass_only,
        'phased_proximal_variants_vcf': args.phased_proximal_variants_vcf,
        'biotypes': args.biotypes,
        'allow_incomplete_transcripts': args.allow_incomplete_transcripts,
        'downstream_sequence_length': downstream_sequence_length,
        'flanking_sequence_length': args.flanking_sequence_length,
        'mutant_only': args.mutant_only,
        'aggregate_report_evaluation': args.aggregate_report_evaluation,
        'input_tsv': args.input_tsv,
        'output_file': args.output_file,
    }
    PvacseqGenerateProteinFasta(**params).execute()


if __name__ == '__main__':
    main()
