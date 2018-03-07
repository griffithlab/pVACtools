import sys
from abc import ABCMeta, abstractmethod
import os
import csv
import datetime
import time

try:
    from .. import lib
except ValueError:
    import lib
from lib.prediction_class import *
from lib.input_file_converter import *
from lib.fasta_generator import *
from lib.output_parser import *
from lib.binding_filter import *
from lib.top_score_filter import *
from lib.filter import *
import shutil
import yaml
import pkg_resources

def status_message(msg):
    print(msg)
    sys.stdout.flush()

class Pipeline(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        self.input_file                  = kwargs['input_file']
        self.input_file_type             = kwargs['input_file_type']
        self.sample_name                 = kwargs['sample_name']
        self.alleles                     = kwargs['alleles']
        self.prediction_algorithms       = kwargs['prediction_algorithms']
        self.output_dir                  = kwargs['output_dir']
        self.iedb_executable             = kwargs.pop('iedb_executable', None)
        self.gene_expn_file              = kwargs.pop('gene_expn_file', None)
        self.transcript_expn_file        = kwargs.pop('transcript_expn_file', None)
        self.normal_snvs_coverage_file   = kwargs.pop('normal_snvs_coverage_file', None)
        self.normal_indels_coverage_file = kwargs.pop('normal_indels_coverage_file', None)
        self.tdna_snvs_coverage_file     = kwargs.pop('tdna_snvs_coverage_file', None)
        self.tdna_indels_coverage_file   = kwargs.pop('tdna_indels_coverage_file', None)
        self.trna_snvs_coverage_file     = kwargs.pop('trna_snvs_coverage_file', None)
        self.trna_indels_coverage_file   = kwargs.pop('trna_indels_coverage_file', None)
        self.net_chop_method             = kwargs.pop('net_chop_method', None)
        self.net_chop_threshold          = kwargs.pop('net_chop_threshold', 0.5)
        self.netmhc_stab                 = kwargs.pop('netmhc_stab', False)
        self.top_result_per_mutation     = kwargs.pop('top_result_per_mutation', False)
        self.top_score_metric            = kwargs.pop('top_score_metric', 'median')
        self.binding_threshold           = kwargs.pop('binding_threshold', 500)
        self.minimum_fold_change         = kwargs.pop('minimum_fold_change', 0)
        self.normal_cov                  = kwargs.pop('normal_cov', None)
        self.normal_vaf                  = kwargs.pop('normal_vaf', None)
        self.tdna_cov                    = kwargs.pop('tdna_cov', None)
        self.tdna_vaf                    = kwargs.pop('tdna_vaf', None)
        self.trna_cov                    = kwargs.pop('trna_cov', None)
        self.trna_vaf                    = kwargs.pop('trna_vaf', None)
        self.expn_val                    = kwargs.pop('expn_val', None)
        self.additional_report_columns   = kwargs.pop('additional_report_columns', None)
        self.fasta_size                  = kwargs.pop('fasta_size', 200)
        self.iedb_retries                = kwargs.pop('iedb_retries', 5)
        self.downstream_sequence_length  = kwargs.pop('downstream_sequence_length', 1000)
        self.keep_tmp_files              = kwargs.pop('keep_tmp_files', False)
        self.exclude_NAs                 = kwargs.pop('exclude_NAs', False)
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
                past_inputs = yaml.load(log_fh)
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

    def converter(self, params):
        converter_types = {
            'vcf'  : 'VcfConverter',
            'bedpe': 'IntegrateConverter',
        }
        converter_type = converter_types[self.input_file_type]
        converter = getattr(sys.modules[__name__], converter_type)
        return converter(**params)

    def fasta_generator(self, params):
        generator_types = {
            'vcf'                   : 'FastaGenerator',
            'bedpe'                 : 'FusionFastaGenerator',
            'pvacvector_input_fasta': 'VectorFastaGenerator',
        }
        generator_type = generator_types[self.input_file_type]
        generator = getattr(sys.modules[__name__], generator_type)
        return generator(**params)

    def output_parser(self, params):
        parser_types = {
            'vcf'  : 'DefaultOutputParser',
            'bedpe': 'FusionOutputParser',
            'pvacvector_input_fasta': 'VectorOutputParser',
        }
        parser_type = parser_types[self.input_file_type]
        parser = getattr(sys.modules[__name__], parser_type)
        return parser(**params)

    def convert_vcf(self):
        status_message("Converting .%s to TSV" % self.input_file_type)
        if os.path.exists(self.tsv_file_path()):
            status_message("TSV file already exists. Skipping.")
            return

        convert_params = {
            'input_file' : self.input_file,
            'output_file': self.tsv_file_path(),
        }
        for attribute in [
            'gene_expn_file',
            'transcript_expn_file',
            'normal_snvs_coverage_file',
            'normal_indels_coverage_file',
            'tdna_snvs_coverage_file',
            'tdna_indels_coverage_file',
            'trna_snvs_coverage_file',
            'trna_indels_coverage_file'
        ]:
            if getattr(self, attribute):
                convert_params[attribute] = getattr(self, attribute)
            else:
                convert_params[attribute] = None

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

    @abstractmethod
    def generate_fasta(self):
        pass

    def split_fasta_basename(self):
        return os.path.join(self.tmp_dir, self.sample_name + "_" + str(self.peptide_sequence_length) + ".fa.split")

    @abstractmethod
    def call_iedb_and_parse_outputs(self, chunks):
        pass

    def combined_parsed_path(self):
        combined_parsed = "%s.combined.parsed.tsv" % self.sample_name
        return os.path.join(self.output_dir, combined_parsed)

    def combined_parsed_outputs(self, split_parsed_output_files):
        status_message("Combining Parsed IEDB Output Files")
        lib.combine_parsed_outputs.main([
            *split_parsed_output_files,
            self.combined_parsed_path(),
            '--top-score-metric', self.top_score_metric,
        ])
        status_message("Completed")

    def binding_filter_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".filtered.binding.tsv")

    def binding_filter(self):
        status_message("Running Binding Filters")
        BindingFilter(
            self.combined_parsed_path(),
            self.binding_filter_out_path(),
            self.binding_threshold,
            self.minimum_fold_change,
            self.top_score_metric,
            self.exclude_NAs,
        ).execute()
        status_message("Completed")

    def coverage_filter_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".filtered.coverage.tsv")

    def coverage_filter(self):
        status_message("Running Coverage Filters")
        if self.input_file_type == 'vcf':
            filter_criteria = []
            filter_criteria.append({'column': "Normal_Depth", 'operator': '>=', 'threshold': self.normal_cov})
            filter_criteria.append({'column': "Normal_VAF", 'operator': '<=', 'threshold': self.normal_vaf})
            filter_criteria.append({'column': "Tumor_DNA_Depth", 'operator': '>=', 'threshold': self.tdna_cov})
            filter_criteria.append({'column': "Tumor_DNA_VAF", 'operator': '>=', 'threshold': self.tdna_vaf})
            filter_criteria.append({'column': "Tumor_RNA_Depth", 'operator': '>=', 'threshold': self.trna_cov})
            filter_criteria.append({'column': "Tumor_RNA_VAF", 'operator': '>=', 'threshold': self.trna_vaf})
            filter_criteria.append({'column': "Gene_Expression", 'operator': '>=', 'threshold': self.expn_val})
            filter_criteria.append({'column': "Transcript_Expression", 'operator': '>=', 'threshold': self.expn_val})
            Filter(self.binding_filter_out_path(), self.coverage_filter_out_path(), filter_criteria, self.exclude_NAs).execute()
        elif self.input_file_type == 'bedpe':
            shutil.copy(self.binding_filter_out_path(), self.coverage_filter_out_path())
        status_message("Completed")

    def top_result_filter_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".filtered.top.tsv")

    def top_result_filter(self):
        status_message("Running Top Score Filter")
        TopScoreFilter(self.coverage_filter_out_path(), self.top_result_filter_out_path(), self.top_score_metric).execute()

    def net_chop_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".chop.tsv")

    def net_chop(self):
        status_message("Submitting remaining epitopes to NetChop")
        lib.net_chop.main([
            self.coverage_filter_out_path(),
            self.net_chop_out_path(),
            '--method',
            self.net_chop_method,
            '--threshold',
            str(self.net_chop_threshold)
        ])
        status_message("Completed")

    def netmhc_stab_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".stab.tsv")

    def call_netmhc_stab(self):
        status_message("Running NetMHCStabPan")
        lib.netmhc_stab.main([
            self.net_chop_out_path(),
            self.netmhc_stab_out_path(),
        ])
        status_message("Completed")

    def final_path(self):
        return os.path.join(self.output_dir, self.sample_name+".final.tsv")

    def execute(self):
        self.print_log()
        self.convert_vcf()

        total_row_count = self.tsv_entry_count()
        if total_row_count == 0:
            if self.input_file_type == 'vcf':
                sys.exit("The TSV file is empty. Please check that the input VCF contains missense, inframe indel, or frameshift mutations.")
            elif self.input_file_type == 'bedpe':
                sys.exit("The TSV file is empty. Please check that the input bedpe file contains fusion entries.")
        chunks = self.split_tsv_file(total_row_count)

        self.generate_fasta(chunks)
        split_parsed_output_files = self.call_iedb_and_parse_outputs(chunks)

        if len(split_parsed_output_files) == 0:
            status_message("No output files were created. Aborting.")
            return

        self.combined_parsed_outputs(split_parsed_output_files)
        self.binding_filter()

        symlinks_to_delete = []
        if (self.gene_expn_file is not None
            or self.transcript_expn_file is not None
            or self.normal_snvs_coverage_file is not None
            or self.normal_indels_coverage_file is not None
            or self.tdna_snvs_coverage_file is not None
            or self.tdna_indels_coverage_file is not None
            or self.trna_snvs_coverage_file is not None
            or self.trna_indels_coverage_file is not None):
            self.coverage_filter()
        else:
            os.symlink(self.binding_filter_out_path(), self.coverage_filter_out_path())
            symlinks_to_delete.append(self.coverage_filter_out_path())

        if self.top_result_per_mutation:
            self.top_result_filter()
        else:
            os.symlink(self.coverage_filter_out_path(), self.top_result_filter_out_path())
            symlinks_to_delete.append(self.top_result_filter_out_path())

        if self.net_chop_method:
            self.net_chop()
        else:
            os.symlink(self.top_result_filter_out_path(), self.net_chop_out_path())
            symlinks_to_delete.append(self.net_chop_out_path())

        if self.netmhc_stab:
            self.call_netmhc_stab()
        else:
            os.symlink(self.net_chop_out_path(), self.netmhc_stab_out_path())
            symlinks_to_delete.append(self.netmhc_stab_out_path())

        shutil.copy(self.netmhc_stab_out_path(), self.final_path())
        for symlink in symlinks_to_delete:
            os.unlink(symlink)


        status_message(
            "\n"
            + "Done: Pipeline finished successfully. File %s contains list of filtered putative neoantigens. " % self.final_path()
        )
        if self.keep_tmp_files is False:
            shutil.rmtree(self.tmp_dir)

class MHCIPipeline(Pipeline):
    def __init__(self, **kwargs):
        Pipeline.__init__(self, **kwargs)
        self.peptide_sequence_length = kwargs.pop('peptide_sequence_length', 21)
        self.epitope_lengths         = kwargs['epitope_lengths']

    def generate_fasta(self, chunks):
        status_message("Generating Variant Peptide FASTA and Key Files")
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            if self.input_file_type == 'pvacvector_input_fasta':
                split_tsv_file_path = self.tsv_file_path()
            else:
                split_tsv_file_path       = "%s_%s" % (self.tsv_file_path(), tsv_chunk)
            split_fasta_file_path     = "%s_%s" % (self.split_fasta_basename(), fasta_chunk)
            if os.path.exists(split_fasta_file_path):
                status_message("Split FASTA file for Entries %s already exists. Skipping." % (fasta_chunk))
                continue
            split_fasta_key_file_path = split_fasta_file_path + '.key'
            status_message("Generating Variant Peptide FASTA and Key Files - Entries %s" % (fasta_chunk))
            generate_fasta_params = {
                'input_file'                : split_tsv_file_path,
                'peptide_sequence_length'   : self.peptide_sequence_length,
                'epitope_length'            : max(self.epitope_lengths),
                'output_file'               : split_fasta_file_path,
                'output_key_file'           : split_fasta_key_file_path,
                'downstream_sequence_length': self.downstream_sequence_length,
            }
            fasta_generator = self.fasta_generator(generate_fasta_params)
            fasta_generator.execute()
        status_message("Completed")

    def call_iedb_and_parse_outputs(self, chunks):
        split_parsed_output_files = []
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            for a in self.alleles:
                for epl in self.epitope_lengths:
                    split_fasta_file_path = "%s_%s"%(self.split_fasta_basename(), fasta_chunk)
                    split_iedb_output_files = []
                    status_message("Processing entries for Allele %s and Epitope Length %s - Entries %s" % (a, epl, fasta_chunk))
                    if os.path.getsize(split_fasta_file_path) == 0:
                        status_message("Fasta file is empty. Skipping")
                        continue
                    for method in self.prediction_algorithms:
                        prediction_class = globals()[method]
                        prediction = prediction_class()
                        iedb_method = prediction.iedb_prediction_method
                        valid_alleles = prediction.valid_allele_names()
                        if a not in valid_alleles:
                            status_message("Allele %s not valid for Method %s. Skipping." % (a, method))
                            continue
                        valid_lengths = prediction.valid_lengths_for_allele(a)
                        if epl not in valid_lengths:
                            status_message("Epitope Length %s is not valid for Method %s and Allele %s. Skipping." % (epl, method, a))
                            continue

                        split_iedb_out = os.path.join(self.tmp_dir, ".".join([self.sample_name, iedb_method, a, str(epl), "tsv_%s" % fasta_chunk]))
                        if os.path.exists(split_iedb_out):
                            status_message("IEDB file for Allele %s and Epitope Length %s with Method %s (Entries %s) already exists. Skipping." % (a, epl, method, fasta_chunk))
                            split_iedb_output_files.append(split_iedb_out)
                            continue
                        status_message("Running IEDB on Allele %s and Epitope Length %s with Method %s - Entries %s" % (a, epl, method, fasta_chunk))

                        if not os.environ.get('TEST_FLAG') or os.environ.get('TEST_FLAG') == '0':
                            if 'last_execute_timestamp' in locals() and not self.iedb_executable:
                                elapsed_time = ( datetime.datetime.now() - last_execute_timestamp ).total_seconds()
                                wait_time = 60 - elapsed_time
                                if wait_time > 0:
                                    time.sleep(wait_time)

                        lib.call_iedb.main([
                            split_fasta_file_path,
                            split_iedb_out,
                            iedb_method,
                            a,
                            '-l', str(epl),
                            '-r', str(self.iedb_retries),
                            '-e', self.iedb_executable,
                        ])
                        last_execute_timestamp = datetime.datetime.now()
                        status_message("Completed")
                        split_iedb_output_files.append(split_iedb_out)

                    split_parsed_file_path = os.path.join(self.tmp_dir, ".".join([self.sample_name, a, str(epl), "parsed", "tsv_%s" % fasta_chunk]))
                    if os.path.exists(split_parsed_file_path):
                        status_message("Parsed Output File for Allele %s and Epitope Length %s (Entries %s) already exists. Skipping" % (a, epl, fasta_chunk))
                        split_parsed_output_files.append(split_parsed_file_path)
                        continue
                    split_fasta_key_file_path = split_fasta_file_path + '.key'

                    if len(split_iedb_output_files) > 0:
                        status_message("Parsing IEDB Output for Allele %s and Epitope Length %s - Entries %s" % (a, epl, fasta_chunk))
                        split_tsv_file_path = "%s_%s" % (self.tsv_file_path(), tsv_chunk)
                        params = {
                            'input_iedb_files'       : split_iedb_output_files,
                            'input_tsv_file'         : split_tsv_file_path,
                            'key_file'               : split_fasta_key_file_path,
                            'output_file'            : split_parsed_file_path,
                        }
                        if self.additional_report_columns and 'sample_name' in self.additional_report_columns:
                            params['sample_name'] = self.sample_name
                        else:
                            params['sample_name'] = None
                        parser = self.output_parser(params)
                        parser.execute()
                        status_message("Completed")
                        split_parsed_output_files.append(split_parsed_file_path)
        return split_parsed_output_files

class MHCIIPipeline(Pipeline):
    def __init__(self, **kwargs):
        Pipeline.__init__(self, **kwargs)
        self.peptide_sequence_length = 31

    def generate_fasta(self, chunks):
        status_message("Generating Variant Peptide FASTA and Key Files")
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            if self.input_file_type == 'pvacvector_input_fasta':
                split_tsv_file_path = self.tsv_file_path()
            else:
                split_tsv_file_path       = "%s_%s" % (self.tsv_file_path(), tsv_chunk)
            split_fasta_file_path     = "%s_%s" % (self.split_fasta_basename(), fasta_chunk)
            if os.path.exists(split_fasta_file_path):
                status_message("Split FASTA file for Entries %s already exists. Skipping." % (fasta_chunk))
                continue
            split_fasta_key_file_path = split_fasta_file_path + '.key'
            status_message("Generating Variant Peptide FASTA and Key Files - Entries %s" % (fasta_chunk))
            generate_fasta_params = {
                'input_file'                : split_tsv_file_path,
                'peptide_sequence_length'   : self.peptide_sequence_length,
                'epitope_length'            : 15,
                'output_file'               : split_fasta_file_path,
                'output_key_file'           : split_fasta_key_file_path,
                'downstream_sequence_length': self.downstream_sequence_length,
            }
            fasta_generator = self.fasta_generator(generate_fasta_params)
            fasta_generator.execute()
        status_message("Completed")

    def call_iedb_and_parse_outputs(self, chunks):
        split_parsed_output_files = []
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            for a in self.alleles:
                split_fasta_file_path = "%s_%s"%(self.split_fasta_basename(), fasta_chunk)
                split_iedb_output_files = []
                status_message("Processing entries for Allele %s - Entries %s" % (a, fasta_chunk))
                if os.path.getsize(split_fasta_file_path) == 0:
                    status_message("Fasta file is empty. Skipping")
                    continue
                for method in self.prediction_algorithms:
                    prediction_class = globals()[method]
                    prediction = prediction_class()
                    iedb_method = prediction.iedb_prediction_method
                    valid_alleles = prediction.valid_allele_names()
                    if a not in valid_alleles:
                        status_message("Allele %s not valid for Method %s. Skipping." % (a, method))
                        continue

                    split_iedb_out = os.path.join(self.tmp_dir, ".".join([self.sample_name, iedb_method, a, "tsv_%s" % fasta_chunk]))
                    if os.path.exists(split_iedb_out):
                        status_message("IEDB file for Allele %s with Method %s (Entries %s) already exists. Skipping." % (a, method, fasta_chunk))
                        split_iedb_output_files.append(split_iedb_out)
                        continue
                    status_message("Running IEDB on Allele %s with Method %s - Entries %s" % (a, method, fasta_chunk))

                    if not os.environ.get('TEST_FLAG') or os.environ.get('TEST_FLAG') == '0':
                        if 'last_execute_timestamp' in locals() and not self.iedb_executable:
                            elapsed_time = ( datetime.datetime.now() - last_execute_timestamp ).total_seconds()
                            wait_time = 60 - elapsed_time
                            if wait_time > 0:
                                time.sleep(wait_time)

                    lib.call_iedb.main([
                        split_fasta_file_path,
                        split_iedb_out,
                        iedb_method,
                        a,
                        '-r', str(self.iedb_retries),
                        '-e', self.iedb_executable,
                    ])
                    last_execute_timestamp = datetime.datetime.now()
                    status_message("Completed")
                    split_iedb_output_files.append(split_iedb_out)

                split_parsed_file_path = os.path.join(self.tmp_dir, ".".join([self.sample_name, a, "parsed", "tsv_%s" % fasta_chunk]))
                if os.path.exists(split_parsed_file_path):
                    status_message("Parsed Output File for Allele %s (Entries %s) already exists. Skipping" % (a, fasta_chunk))
                    split_parsed_output_files.append(split_parsed_file_path)
                    continue
                split_fasta_key_file_path = split_fasta_file_path + '.key'

                if len(split_iedb_output_files) > 0:
                    status_message("Parsing IEDB Output for Allele %s - Entries %s" % (a, fasta_chunk))
                    split_tsv_file_path = "%s_%s" % (self.tsv_file_path(), tsv_chunk)
                    params = {
                        'input_iedb_files'       : split_iedb_output_files,
                        'input_tsv_file'         : split_tsv_file_path,
                        'key_file'               : split_fasta_key_file_path,
                        'output_file'            : split_parsed_file_path,
                    }
                    if self.additional_report_columns and 'sample_name' in self.additional_report_columns:
                        params['sample_name'] = self.sample_name
                    else:
                        params['sample_name'] = None
                    parser = self.output_parser(params)
                    parser.execute()
                    status_message("Completed")
                    split_parsed_output_files.append(split_parsed_file_path)

        return split_parsed_output_files
