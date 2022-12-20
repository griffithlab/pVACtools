import sys
from abc import ABCMeta, abstractmethod
import os
import csv
import datetime
import time
import shutil
import copy
import yaml
import pkg_resources
import pymp
from threading import Lock
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from collections import OrderedDict
import logging

from pvactools.lib.prediction_class import *
from pvactools.lib.input_file_converter import VcfConverter
from pvactools.lib.fasta_generator import FastaGenerator, VectorFastaGenerator
from pvactools.lib.output_parser import DefaultOutputParser, UnmatchedSequencesOutputParser
from pvactools.lib.post_processor import PostProcessor
import pvactools.lib.call_iedb
import pvactools.lib.combine_parsed_outputs

def status_message(msg):
    print(msg)
    sys.stdout.flush()

class Pipeline(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        for (k,v) in kwargs.items():
           setattr(self, k, v)
        self.proximal_variants_file      = None
        tmp_dir = os.path.join(self.output_dir, 'tmp')
        os.makedirs(tmp_dir, exist_ok=True)
        self.tmp_dir = tmp_dir

    def log_dir(self):
        dir = os.path.join(self.output_dir, 'log')
        os.makedirs(dir, exist_ok=True)
        return dir

    def print_log(self):
        log_file = os.path.join(self.log_dir(), 'inputs.yml')
        if os.path.exists(log_file):
            with open(log_file, 'r') as log_fh:
                past_inputs = yaml.load(log_fh, Loader=yaml.FullLoader)
                current_inputs = self.__dict__
                current_inputs['pvactools_version'] = pkg_resources.get_distribution("pvactools").version
                if past_inputs['pvactools_version'] != current_inputs['pvactools_version']:
                    status_message(
                        "Restart to be executed with a different pVACtools version:\n" +
                        "Past version: %s\n" % past_inputs['pvactools_version'] +
                        "Current version: %s" % current_inputs['pvactools_version']
                    )
                for key in current_inputs.keys():
                    if key == 'pvactools_version' or key == 'pvacseq_version':
                        continue
                    if key not in past_inputs.keys() and current_inputs[key] is not None:
                        sys.exit(
                            "Restart inputs are different from past inputs: \n" +
                            "Additional input: %s - %s\n" % (key, current_inputs[key]) +
                            "Aborting."
                        )
                    elif current_inputs[key] != past_inputs[key]:
                        sys.exit(
                            "Restart inputs are different from past inputs: \n" +
                            "Past input: %s - %s\n" % (key, past_inputs[key]) +
                            "Current input: %s - %s\n" % (key, current_inputs[key]) +
                            "Aborting."
                        )
        else:
            with open(log_file, 'w') as log_fh:
                inputs = self.__dict__
                inputs['pvactools_version'] = pkg_resources.get_distribution("pvactools").version
                yaml.dump(inputs, log_fh, default_flow_style=False)

    def tsv_file_path(self):
        if self.input_file_type == 'pvacvector_input_fasta':
            return self.input_file
        else:
            tsv_file = self.sample_name + '.tsv'
            return os.path.join(self.output_dir, tsv_file)

    def fasta_file_path(self):
        fasta_file = self.sample_name + '.fasta'
        return os.path.join(self.output_dir, fasta_file)

    def net_chop_fasta_file_path(self):
        net_chop_fasta_file = self.sample_name + '.net_chop.fa'
        return os.path.join(self.output_dir, net_chop_fasta_file)

    def converter(self, params):
        converter_types = {
            'vcf'  : 'VcfConverter',
        }
        converter_type = converter_types[self.input_file_type]
        converter = getattr(sys.modules[__name__], converter_type)
        return converter(**params)

    def fasta_generator(self, params):
        generator_types = {
            'vcf'                   : 'FastaGenerator',
            'pvacvector_input_fasta': 'VectorFastaGenerator',
        }
        generator_type = generator_types[self.input_file_type]
        generator = getattr(sys.modules[__name__], generator_type)
        return generator(**params)

    def output_parser(self, params):
        parser_types = {
            'vcf'  : 'DefaultOutputParser',
            'pvacvector_input_fasta': 'UnmatchedSequencesOutputParser',
            'fasta': 'UnmatchedSequencesOutputParser',
        }
        parser_type = parser_types[self.input_file_type]
        parser = getattr(sys.modules[__name__], parser_type)
        return parser(**params)

    def generate_combined_fasta(self, fasta_path, epitope_flank_length=0):
        params = [
            self.input_file,
            str(epitope_flank_length + max(self.epitope_lengths) - 1),
            fasta_path,
        ]
        import pvactools.tools.pvacseq.generate_protein_fasta as generate_combined_fasta
        params.extend(["--sample-name", self.sample_name])
        if self.phased_proximal_variants_vcf is not None:
            params.extend(["--phased-proximal-variants-vcf", self.phased_proximal_variants_vcf])
        if self.downstream_sequence_length is not None:
            params.extend(["-d", str(self.downstream_sequence_length)])
        else:
            params.extend(["-d", 'full'])
        if self.pass_only:
            params.extend(["--pass-only"])
        generate_combined_fasta.main(params)
        os.unlink("{}.manufacturability.tsv".format(fasta_path))

    def convert_vcf(self):
        status_message("Converting .%s to TSV" % self.input_file_type)
        if os.path.exists(self.tsv_file_path()):
            status_message("TSV file already exists. Skipping.")
            return

        convert_params = {
            'input_file' : self.input_file,
            'output_file': self.tsv_file_path(),
            'sample_name': self.sample_name,
            'pass_only': self.pass_only,
        }
        if self.normal_sample_name is not None:
            convert_params['normal_sample_name'] = self.normal_sample_name
        if self.phased_proximal_variants_vcf is not None:
            convert_params['proximal_variants_vcf'] = self.phased_proximal_variants_vcf
            proximal_variants_tsv = os.path.join(self.output_dir, self.sample_name + '.proximal_variants.tsv')
            convert_params['proximal_variants_tsv'] = proximal_variants_tsv
            self.proximal_variants_file = proximal_variants_tsv
            convert_params['flanking_bases'] = max(self.epitope_lengths) * 4

        converter = self.converter(convert_params)
        converter.execute()
        print("Completed")

    def tsv_entry_count(self):
        with open(self.tsv_file_path()) as tsv_file:
            reader  = csv.DictReader(tsv_file, delimiter='\t')
            row_count = 0
            for row in reader:
                row_count += 1
        return row_count

    def split_tsv_file(self, total_row_count):
        status_message("Splitting TSV into smaller chunks")
        tsv_size = self.fasta_size / 2
        chunks = []
        with open(self.tsv_file_path(), 'r') as tsv_file:
            reader      = csv.DictReader(tsv_file, delimiter='\t')
            row_count   = 1
            split_start = row_count
            split_end   = split_start + tsv_size - 1
            if split_end > total_row_count:
                split_end = total_row_count
            status_message("Splitting TSV into smaller chunks - Entries %d-%d" % (split_start, split_end))
            split_tsv_file_path = "%s_%d-%d" % (self.tsv_file_path(), split_start, split_end)
            chunks.append([split_start, split_end])
            if os.path.exists(split_tsv_file_path):
                status_message("Split TSV file for Entries %d-%d already exists. Skipping." % (split_start, split_end))
                skip = 1
            else:
                split_tsv_file      = open(split_tsv_file_path, 'w')
                split_tsv_writer    = csv.DictWriter(split_tsv_file, delimiter='\t', fieldnames = reader.fieldnames)
                split_tsv_writer.writeheader()
                skip = 0
            for row in reader:
                if skip == 0:
                    split_tsv_writer.writerow(row)
                if row_count == total_row_count:
                    break
                if row_count % tsv_size == 0:
                    if skip == 0:
                        split_tsv_file.close()
                    split_start = row_count + 1
                    split_end   = split_start + tsv_size - 1
                    if split_end > total_row_count:
                        split_end = total_row_count
                    status_message("Splitting TSV into smaller chunks - Entries %d-%d" % (split_start, split_end))
                    split_tsv_file_path = "%s_%d-%d" % (self.tsv_file_path(), split_start, split_end)
                    chunks.append([split_start, split_end])
                    if os.path.exists(split_tsv_file_path):
                        status_message("Split TSV file for Entries %d-%d already exists. Skipping." % (split_start, split_end))
                        skip = 1
                    else:
                        split_tsv_file      = open(split_tsv_file_path, 'w')
                        split_tsv_writer    = csv.DictWriter(split_tsv_file, delimiter='\t', fieldnames = reader.fieldnames)
                        split_tsv_writer.writeheader()
                        skip = 0
                row_count += 1
            if skip == 0:
                split_tsv_file.close()
        status_message("Completed")
        return chunks

    def generate_fasta(self, chunks):
        status_message("Generating Variant Peptide FASTA and Key Files")
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            generate_fasta_params = {
                'downstream_sequence_length': self.downstream_sequence_length,
                'proximal_variants_file'    : self.proximal_variants_file,
            }
            if self.input_file_type == 'pvacvector_input_fasta':
                split_fasta_file_path = "{}_{}".format(self.split_fasta_basename(None), fasta_chunk)
                generate_fasta_params['input_file'] = self.tsv_file_path()
                generate_fasta_params['output_file_prefix'] = split_fasta_file_path
                generate_fasta_params['epitope_lengths'] = self.epitope_lengths
                generate_fasta_params['spacers'] = self.spacers
                status_message("Generating Variant Peptide FASTA and Key Files - Entries %s" % (fasta_chunk))
                fasta_generator = self.fasta_generator(generate_fasta_params)
                fasta_generator.execute()
            else:
                for epitope_length in self.epitope_lengths:
                    split_fasta_file_path = "{}_{}".format(self.split_fasta_basename(epitope_length), fasta_chunk)
                    if os.path.exists(split_fasta_file_path):
                        status_message("Split FASTA file for Epitope Length {} - Entries {} already exists. Skipping.".format(epitope_length, fasta_chunk))
                        continue
                    split_fasta_key_file_path = split_fasta_file_path + '.key'
                    generate_fasta_params['input_file'] = "%s_%s" % (self.tsv_file_path(), tsv_chunk)
                    generate_fasta_params['epitope_length'] = epitope_length
                    generate_fasta_params['flanking_sequence_length'] = epitope_length - 1
                    generate_fasta_params['output_file'] = split_fasta_file_path
                    generate_fasta_params['output_key_file'] = split_fasta_key_file_path
                    status_message("Generating Variant Peptide FASTA and Key Files - Epitope Length {} - Entries {}".format(epitope_length, fasta_chunk))
                    fasta_generator = self.fasta_generator(generate_fasta_params)
                    fasta_generator.execute()
        status_message("Completed")

    def split_fasta_basename(self, epitope_length):
        if epitope_length is None:
            return os.path.join(self.tmp_dir, "{}.fa.split".format(self.sample_name))
        else:
            return os.path.join(self.tmp_dir, "{}.{}.fa.split".format(self.sample_name, epitope_length))

    def call_iedb(self, chunks):
        alleles = self.alleles
        epitope_lengths = self.epitope_lengths
        prediction_algorithms = self.prediction_algorithms
        argument_sets = []
        warning_messages = []
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            if self.input_file_type == 'fasta':
                fasta_chunk = tsv_chunk
            else:
                fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            for a in alleles:
                for epl in epitope_lengths:
                    if self.input_file_type == 'pvacvector_input_fasta':
                        split_fasta_file_path = "{}_1-2.{}.tsv".format(self.split_fasta_basename(None), epl)
                    else:
                        split_fasta_file_path = "%s_%s"%(self.split_fasta_basename(epl), fasta_chunk)
                    if os.path.getsize(split_fasta_file_path) == 0:
                        msg = "Fasta file {} is empty. Skipping".format(split_fasta_file_path)
                        if msg not in warning_messages:
                            warning_messages.append(msg)
                        continue
                    #begin of per-algorithm processing
                    for method in prediction_algorithms:
                        prediction_class = globals()[method]
                        prediction = prediction_class()
                        if hasattr(prediction, 'iedb_prediction_method'):
                            iedb_method = prediction.iedb_prediction_method
                        else:
                            iedb_method = method
                        valid_alleles = prediction.valid_allele_names()
                        if a not in valid_alleles:
                            msg = "Allele %s not valid for Method %s. Skipping." % (a, method)
                            if msg not in warning_messages:
                                warning_messages.append(msg)
                            continue
                        valid_lengths = prediction.valid_lengths_for_allele(a)
                        if epl not in valid_lengths:
                            msg = "Epitope Length %s is not valid for Method %s and Allele %s. Skipping." % (epl, method, a)
                            if msg not in warning_messages:
                                warning_messages.append(msg)
                            continue

                        split_iedb_out = os.path.join(self.tmp_dir, ".".join([self.sample_name, iedb_method, a, str(epl), "tsv_%s" % fasta_chunk]))
                        if os.path.exists(split_iedb_out):
                            msg = "Prediction file for Allele %s and Epitope Length %s with Method %s (Entries %s) already exists. Skipping." % (a, epl, method, fasta_chunk)
                            if msg not in warning_messages:
                                warning_messages.append(msg)
                            continue
                        arguments = [
                            split_fasta_file_path,
                            split_iedb_out,
                            method,
                            a,
                            '-r', str(self.iedb_retries),
                            '-e', self.iedb_executable,
                            '-l', str(epl),
                        ]
                        argument_sets.append(arguments)

        for msg in warning_messages:
            status_message(msg)

        with pymp.Parallel(self.n_threads) as p:
            for index in p.range(len(argument_sets)):
                arguments = argument_sets[index]
                a = arguments[3]
                method = arguments[2]
                filename = arguments[1]
                epl = arguments[9]
                p.print("Making binding predictions on Allele %s and Epitope Length %s with Method %s - File %s" % (a, epl, method, filename))
                pvactools.lib.call_iedb.main(arguments)
                p.print("Making binding predictions on Allele %s and Epitope Length %s with Method %s - File %s - Completed" % (a, epl, method, filename))

    def parse_outputs(self, chunks):
        split_parsed_output_files = []
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            if self.input_file_type == 'fasta':
                fasta_chunk = tsv_chunk
            else:
                fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            for a in self.alleles:
                for epl in self.epitope_lengths:
                    split_iedb_output_files = []
                    status_message("Parsing binding predictions for Allele %s and Epitope Length %s - Entries %s" % (a, epl, fasta_chunk))
                    for method in self.prediction_algorithms:
                        prediction_class = globals()[method]
                        prediction = prediction_class()
                        if hasattr(prediction, 'iedb_prediction_method'):
                            iedb_method = prediction.iedb_prediction_method
                        else:
                            iedb_method = method
                        valid_alleles = prediction.valid_allele_names()
                        if a not in valid_alleles:
                            continue
                        valid_lengths = prediction.valid_lengths_for_allele(a)
                        if epl not in valid_lengths:
                            continue
                        split_iedb_out = os.path.join(self.tmp_dir, ".".join([self.sample_name, iedb_method, a, str(epl), "tsv_%s" % fasta_chunk]))
                        if os.path.exists(split_iedb_out):
                            split_iedb_output_files.append(split_iedb_out)

                    split_parsed_file_path = os.path.join(self.tmp_dir, ".".join([self.sample_name, a, str(epl), "parsed", "tsv_%s" % fasta_chunk]))
                    if os.path.exists(split_parsed_file_path):
                        status_message("Parsed Output File for Allele %s and Epitope Length %s (Entries %s) already exists. Skipping" % (a, epl, fasta_chunk))
                        split_parsed_output_files.append(split_parsed_file_path)
                        continue
                    if self.input_file_type == 'pvacvector_input_fasta':
                        split_fasta_file_path = "{}_1-2.{}.tsv".format(self.split_fasta_basename(None), epl)
                    else:
                        split_fasta_file_path = "%s_%s"%(self.split_fasta_basename(epl), fasta_chunk)
                    split_fasta_key_file_path = split_fasta_file_path + '.key'

                    if len(split_iedb_output_files) > 0:
                        status_message("Parsing prediction file for Allele %s and Epitope Length %s - Entries %s" % (a, epl, fasta_chunk))
                        split_tsv_file_path = "%s_%s" % (self.tsv_file_path(), tsv_chunk)
                        params = {
                            'input_iedb_files'       : split_iedb_output_files,
                            'input_tsv_file'         : split_tsv_file_path,
                            'key_file'               : split_fasta_key_file_path,
                            'output_file'            : split_parsed_file_path,
                        }
                        params['sample_name'] = self.sample_name
                        if self.additional_report_columns and 'sample_name' in self.additional_report_columns:
                            params['add_sample_name_column'] = True 
                        parser = self.output_parser(params)
                        parser.execute()
                        status_message("Parsing prediction file for Allele %s and Epitope Length %s - Entries %s - Completed" % (a, epl, fasta_chunk))

                        split_parsed_output_files.append(split_parsed_file_path)
        return split_parsed_output_files

    def combined_parsed_path(self):
        combined_parsed = "%s.all_epitopes.tsv" % self.sample_name
        return os.path.join(self.output_dir, combined_parsed)

    def combined_parsed_outputs(self, split_parsed_output_files):
        status_message("Combining Parsed Prediction Files")
        params = [
            *split_parsed_output_files,
            self.combined_parsed_path(),
            '--top-score-metric', self.top_score_metric,
        ]
        if self.input_file_type == 'fasta':
            params.extend(['--file-type', 'pVACbind'])
        pvactools.lib.combine_parsed_outputs.main(params)
        status_message("Completed")

    def final_path(self):
        return os.path.join(self.output_dir, self.sample_name+".filtered.tsv")

    def execute(self):
        self.print_log()
        self.convert_vcf()
        if self.input_file_type != 'pvacvector_input_fasta':
            self.generate_combined_fasta(self.fasta_file_path())

        total_row_count = self.tsv_entry_count()
        if total_row_count == 0:
            print("The TSV file is empty. Please check that the input VCF contains missense, inframe indel, or frameshift mutations.")
            return
        chunks = self.split_tsv_file(total_row_count)

        self.generate_fasta(chunks)
        self.call_iedb(chunks)
        split_parsed_output_files = self.parse_outputs(chunks)

        if len(split_parsed_output_files) == 0:
            status_message("No output files were created. Aborting.")
            return

        self.combined_parsed_outputs(split_parsed_output_files)

        post_processing_params = copy.copy(vars(self))
        post_processing_params['input_file'] = self.combined_parsed_path()
        post_processing_params['file_type'] = self.input_file_type
        post_processing_params['filtered_report_file'] = self.final_path()
        post_processing_params['fasta'] = self.fasta_file_path()
        post_processing_params['run_manufacturability_metrics'] = True
        if self.input_file_type == 'vcf':
            post_processing_params['file_type'] = 'pVACseq'
            post_processing_params['run_coverage_filter'] = True
            post_processing_params['run_transcript_support_level_filter'] = True
        else:
            post_processing_params['run_coverage_filter'] = False
            post_processing_params['run_transcript_support_level_filter'] = False
        if self.net_chop_method:
            net_chop_epitope_flank_length = 9
            self.generate_combined_fasta(self.net_chop_fasta_file_path(), net_chop_epitope_flank_length)
            post_processing_params['net_chop_fasta'] = self.net_chop_fasta_file_path()
            post_processing_params['run_net_chop'] = True
        else:
            post_processing_params['run_net_chop'] = False
        if self.netmhc_stab:
            post_processing_params['run_netmhc_stab'] = True
        else:
            post_processing_params['run_netmhc_stab'] = False
        PostProcessor(**post_processing_params).execute()

        if self.keep_tmp_files is False:
            shutil.rmtree(self.tmp_dir, ignore_errors=True)

class PvacbindPipeline(Pipeline):
    def fasta_entry_count(self):
        with open(self.input_file) as fasta_file:
            row_count = 0
            for row in fasta_file:
                if not row.startswith('>'):
                    row_count += 1
        return row_count

    def fasta_basename(self, length):
        return os.path.join(self.tmp_dir, "{}.{}.fa".format(self.sample_name, length))

    def split_fasta_basename(self, length):
        return "{}.split".format(self.fasta_basename(length))

    def uniquify_records(self, records):
        fasta_sequences = OrderedDict()
        ids = []
        for record in records:
            if record.id in ids:
                raise Exception("Duplicate fasta header {}. Please ensure that the input FASTA uses unique headers.".format(record.id))
            fasta_sequences.setdefault(str(record.seq), []).append(record.id)
            ids.append(record.id)
        count = 1
        uniq_records = []
        keys = {}
        for sequence, ids in fasta_sequences.items():
            record = SeqRecord(Seq(sequence, IUPAC.protein), id=str(count), description=str(count))
            uniq_records.append(record)
            keys[count] = ids
            count += 1
        return (uniq_records, keys)

    def create_per_length_fasta_and_process_stops(self, length):
        stop_chars = set('X*')
        supported_amino_acids = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        records = []
        for record in SeqIO.parse(self.input_file, "fasta"):
            sequence = str(record.seq).upper()
            x_index = sequence.index('X') if 'X' in sequence else len(sequence)
            star_index = sequence.index('*') if '*' in sequence else len(sequence)
            sequence = sequence[0:min(x_index, star_index)]
            if not all([c in supported_amino_acids for c in sequence]):
                logging.warning("Record {} contains unsupported amino acids. Skipping.".format(record.id))
                continue
            if len(sequence) >= length:
                record.seq = Seq(sequence, IUPAC.protein)
                records.append(record)
        SeqIO.write(records, self.fasta_basename(length), "fasta")

    def split_fasta_file(self, length):
        fasta_entry_count = self.fasta_entry_count()
        status_message("Splitting FASTA into smaller chunks")
        chunks = []
        peptides = []
        row_count   = 1
        split_start = row_count
        split_end   = split_start + self.fasta_size - 1
        if split_end > fasta_entry_count:
            split_end = fasta_entry_count
        status_message("Splitting FASTA into smaller chunks - Entries %d-%d" % (split_start, split_end))
        split_fasta_file_path = "%s_%d-%d" % (self.split_fasta_basename(length), split_start, split_end)
        split_fasta_key_file_path = "{}.key".format(split_fasta_file_path)
        chunks.append([split_start, split_end])
        if os.path.exists(split_fasta_file_path):
            status_message("Split FASTA file for Entries %d-%d already exists. Skipping." % (split_start, split_end))
            skip = 1
        else:
            split_fasta_records = []
            skip = 0
        for record in SeqIO.parse(self.fasta_basename(length), "fasta"):
            if skip == 0:
                split_fasta_records.append(record)
            if row_count == fasta_entry_count:
                break
            if row_count % self.fasta_size == 0:
                if skip == 0:
                    (uniq_records, keys) = self.uniquify_records(split_fasta_records)
                    with open(split_fasta_file_path, 'w') as split_fasta_file:
                        SeqIO.write(uniq_records, split_fasta_file, "fasta")
                    with open(split_fasta_key_file_path, 'w') as split_fasta_key_file:
                        yaml.dump(keys, split_fasta_key_file, default_flow_style=False)
                    split_fasta_file.close()
                split_start = row_count + 1
                split_end   = split_start + self.fasta_size - 1
                if split_end > fasta_entry_count:
                    split_end = fasta_entry_count
                status_message("Splitting FASTA into smaller chunks - Entries %d-%d" % (split_start, split_end))
                split_fasta_file_path = "%s_%d-%d" % (self.split_fasta_basename(length), split_start, split_end)
                split_fasta_key_file_path = "{}.key".format(split_fasta_file_path)
                chunks.append([split_start, split_end])
                if os.path.exists(split_fasta_file_path):
                    status_message("Split FASTA file for Entries %d-%d already exists. Skipping." % (split_start, split_end))
                    skip = 1
                else:
                    split_fasta_records = []
                    skip = 0
            row_count += 1
        if skip == 0:
            (uniq_records, keys) = self.uniquify_records(split_fasta_records)
            with open(split_fasta_file_path, 'w') as split_fasta_file:
                SeqIO.write(uniq_records, split_fasta_file, "fasta")
            with open(split_fasta_key_file_path, 'w') as split_fasta_key_file:
                yaml.dump(keys, split_fasta_key_file, default_flow_style=False)
        status_message("Completed")
        return chunks

    def call_iedb(self, chunks, length):
        alleles = self.alleles
        prediction_algorithms = self.prediction_algorithms
        argument_sets = []
        warning_messages = []
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            if self.input_file_type == 'fasta':
                fasta_chunk = tsv_chunk
            else:
                fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            for a in alleles:
                split_fasta_file_path = "%s_%s"%(self.split_fasta_basename(length), fasta_chunk)
                if os.path.getsize(split_fasta_file_path) == 0:
                    msg = "Fasta file {} is empty. Skipping".format(split_fasta_file_path)
                    if msg not in warning_messages:
                        warning_messages.append(msg)
                    continue
                #begin of per-algorithm processing
                for method in prediction_algorithms:
                    prediction_class = globals()[method]
                    prediction = prediction_class()
                    if hasattr(prediction, 'iedb_prediction_method'):
                        iedb_method = prediction.iedb_prediction_method
                    else:
                        iedb_method = method
                    valid_alleles = prediction.valid_allele_names()
                    if a not in valid_alleles:
                        msg = "Allele %s not valid for Method %s. Skipping." % (a, method)
                        if msg not in warning_messages:
                            warning_messages.append(msg)
                        continue
                    valid_lengths = prediction.valid_lengths_for_allele(a)
                    if length not in valid_lengths:
                        msg = "Epitope Length %s is not valid for Method %s and Allele %s. Skipping." % (length, method, a)
                        if msg not in warning_messages:
                            warning_messages.append(msg)
                        continue

                    split_iedb_out = os.path.join(self.tmp_dir, ".".join([self.sample_name, iedb_method, a, str(length), "tsv_%s" % fasta_chunk]))
                    if os.path.exists(split_iedb_out):
                        msg = "Prediction file for Allele %s and Epitope Length %s with Method %s (Entries %s) already exists. Skipping." % (a, length, method, fasta_chunk)
                        if msg not in warning_messages:
                            warning_messages.append(msg)
                        continue
                    arguments = [
                        split_fasta_file_path,
                        split_iedb_out,
                        method,
                        a,
                        '-r', str(self.iedb_retries),
                        '-e', self.iedb_executable,
                        '-l', str(length),
                    ]
                    argument_sets.append(arguments)

        for msg in warning_messages:
            status_message(msg)

        with pymp.Parallel(self.n_threads) as p:
            for index in p.range(len(argument_sets)):
                arguments = argument_sets[index]
                a = arguments[3]
                method = arguments[2]
                filename = arguments[1]
                epl = arguments[9]
                p.print("Making binding predictions on Allele %s and Epitope Length %s with Method %s - File %s" % (a, epl, method, filename))
                pvactools.lib.call_iedb.main(arguments)
                p.print("Making binding predictions on Allele %s and Epitope Length %s with Method %s - File %s - Completed" % (a, epl, method, filename))

    def parse_outputs(self, chunks, length):
        split_parsed_output_files = []
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            if self.input_file_type == 'fasta':
                fasta_chunk = tsv_chunk
            else:
                fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            for a in self.alleles:
                split_iedb_output_files = []
                status_message("Parsing binding predictions for Allele %s and Epitope Length %s - Entries %s" % (a, length, fasta_chunk))
                for method in self.prediction_algorithms:
                    prediction_class = globals()[method]
                    prediction = prediction_class()
                    if hasattr(prediction, 'iedb_prediction_method'):
                        iedb_method = prediction.iedb_prediction_method
                    else:
                        iedb_method = method
                    valid_alleles = prediction.valid_allele_names()
                    if a not in valid_alleles:
                        continue
                    valid_lengths = prediction.valid_lengths_for_allele(a)
                    if length not in valid_lengths:
                        continue
                    split_iedb_out = os.path.join(self.tmp_dir, ".".join([self.sample_name, iedb_method, a, str(length), "tsv_%s" % fasta_chunk]))
                    if os.path.exists(split_iedb_out):
                        split_iedb_output_files.append(split_iedb_out)

                split_parsed_file_path = os.path.join(self.tmp_dir, ".".join([self.sample_name, a, str(length), "parsed", "tsv_%s" % fasta_chunk]))
                if os.path.exists(split_parsed_file_path):
                    status_message("Parsed Output File for Allele %s and Epitope Length %s (Entries %s) already exists. Skipping" % (a, length, fasta_chunk))
                    split_parsed_output_files.append(split_parsed_file_path)
                    continue
                split_fasta_file_path = "%s_%s"%(self.split_fasta_basename(length), fasta_chunk)
                split_fasta_key_file_path = split_fasta_file_path + '.key'

                if len(split_iedb_output_files) > 0:
                    status_message("Parsing prediction file for Allele %s and Epitope Length %s - Entries %s" % (a, length, fasta_chunk))
                    split_tsv_file_path = "%s_%s" % (self.tsv_file_path(), tsv_chunk)
                    params = {
                        'input_iedb_files'       : split_iedb_output_files,
                        'input_tsv_file'         : split_tsv_file_path,
                        'key_file'               : split_fasta_key_file_path,
                        'output_file'            : split_parsed_file_path,
                    }
                    params['sample_name'] = self.sample_name
                    if self.additional_report_columns and 'sample_name' in self.additional_report_columns:
                        params['add_sample_name_column'] = True 
                    parser = self.output_parser(params)
                    parser.execute()
                    status_message("Parsing prediction file for Allele %s and Epitope Length %s - Entries %s - Completed" % (a, length, fasta_chunk))

                    split_parsed_output_files.append(split_parsed_file_path)
        return split_parsed_output_files

    def execute(self):
        self.print_log()

        split_parsed_output_files = []
        for length in self.epitope_lengths:
            self.create_per_length_fasta_and_process_stops(length)
            chunks = self.split_fasta_file(length)
            self.call_iedb(chunks, length)
            split_parsed_output_files.extend(self.parse_outputs(chunks, length))

        if len(split_parsed_output_files) == 0:
            status_message("No output files were created. Aborting.")
            return

        self.combined_parsed_outputs(split_parsed_output_files)

        if not self.run_post_processor:
            return

        post_processing_params = copy.copy(vars(self))
        post_processing_params['input_file'] = self.combined_parsed_path()
        post_processing_params['file_type'] = 'pVACbind'
        post_processing_params['filtered_report_file'] = self.final_path()
        post_processing_params['run_coverage_filter'] = False
        post_processing_params['run_transcript_support_level_filter'] = False
        post_processing_params['minimum_fold_change'] = None
        post_processing_params['run_manufacturability_metrics'] = True
        post_processing_params['fasta'] = self.input_file
        if self.net_chop_method:
            post_processing_params['net_chop_fasta'] = self.net_chop_fasta
            post_processing_params['run_net_chop'] = True
        else:
            post_processing_params['run_net_chop'] = False
        if self.netmhc_stab:
            post_processing_params['run_netmhc_stab'] = True
        else:
            post_processing_params['run_netmhc_stab'] = False
        PostProcessor(**post_processing_params).execute()

        if self.keep_tmp_files is False:
            shutil.rmtree(self.tmp_dir, ignore_errors=True)
