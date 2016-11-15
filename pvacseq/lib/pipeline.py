import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

from abc import ABCMeta, abstractmethod
import os
import csv

try:
    from .. import lib
except ValueError:
    import lib
from lib.prediction_class import *
import shutil

def status_message(msg):
    print(msg)
    sys.stdout.flush()

class Pipeline(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        self.input_file                  = kwargs['input_file']
        self.sample_name                 = kwargs['sample_name']
        self.alleles                     = kwargs['alleles']
        self.prediction_algorithms       = kwargs['prediction_algorithms']
        self.output_dir                  = kwargs['output_dir']
        self.iedb_executable             = kwargs['iedb_executable']
        self.gene_expn_file              = kwargs['gene_expn_file']
        self.transcript_expn_file        = kwargs['transcript_expn_file']
        self.normal_snvs_coverage_file   = kwargs['normal_snvs_coverage_file']
        self.normal_indels_coverage_file = kwargs['normal_indels_coverage_file']
        self.tdna_snvs_coverage_file     = kwargs['tdna_snvs_coverage_file']
        self.tdna_indels_coverage_file   = kwargs['tdna_indels_coverage_file']
        self.trna_snvs_coverage_file     = kwargs['trna_snvs_coverage_file']
        self.trna_indels_coverage_file   = kwargs['trna_indels_coverage_file']
        self.net_chop_method             = kwargs['net_chop_method']
        self.net_chop_threshold          = kwargs['net_chop_threshold']
        self.netmhc_stab                 = kwargs['netmhc_stab']
        self.top_result_per_mutation     = kwargs['top_result_per_mutation']
        self.top_score_metric            = kwargs['top_score_metric']
        self.binding_threshold           = kwargs['binding_threshold']
        self.minimum_fold_change         = kwargs['minimum_fold_change']
        self.normal_cov                  = kwargs['normal_cov']
        self.normal_vaf                  = kwargs['normal_vaf']
        self.tdna_cov                    = kwargs['tdna_cov']
        self.tdna_vaf                    = kwargs['tdna_vaf']
        self.trna_cov                    = kwargs['trna_cov']
        self.trna_vaf                    = kwargs['trna_vaf']
        self.expn_val                    = kwargs['expn_val']
        self.fasta_size                  = kwargs['fasta_size']
        self.iedb_retries                = kwargs['iedb_retries']
        self.downstream_sequence_length  = kwargs['downstream_sequence_length']
        self.keep_tmp_files              = kwargs['keep_tmp_files']
        tmp_dir = os.path.join(self.output_dir, 'tmp')
        os.makedirs(tmp_dir, exist_ok=True)
        self.tmp_dir = tmp_dir

    def tsv_file_path(self):
        tsv_file = self.sample_name + '.tsv'
        return os.path.join(self.output_dir, tsv_file)

    def convert_vcf(self):
        status_message("Converting VCF to TSV")
        if os.path.exists(self.tsv_file_path()):
            status_message("TSV file already exists. Skipping.")
            return

        convert_params = [
            self.input_file,
            self.tsv_file_path(),
        ]
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
                param = '--' + attribute
                param = param.replace('_', '-')
                convert_params.extend([param, getattr(self, attribute)])

        lib.convert_vcf.main(convert_params)
        status_message("Completed")

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
            self.combined_parsed_path()
        ])
        status_message("Completed")

    def binding_filter_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".filtered.binding.tsv")

    def binding_filter(self):
        status_message("Running Binding Filters")
        lib.binding_filter.main(
            [
                self.combined_parsed_path(),
                self.binding_filter_out_path(),
                '-c', str(self.minimum_fold_change),
                '-b', str(self.binding_threshold),
                '-m', str(self.top_score_metric),
            ]
        )
        status_message("Completed")

    def coverage_filter_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".filtered.coverage.tsv")

    def coverage_filter(self):
        status_message("Running Coverage Filters")
        coverage_params = [
            self.binding_filter_out_path(),
            self.coverage_filter_out_path(),
            '--expn-val', str(self.expn_val),
        ]
        for attribute in [
            'expn_val',
            'normal_cov',
            'normal_vaf',
            'tdna_cov',
            'tdna_vaf',
            'trna_cov',
            'trna_vaf',
        ]:
            if getattr(self, attribute):
                param = '--' + attribute
                param = param.replace('_', '-')
                coverage_params.extend([param, str(getattr(self, attribute))])
        lib.coverage_filter.main(coverage_params)
        status_message("Completed")

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
        self.convert_vcf()

        total_row_count = self.tsv_entry_count()
        if total_row_count == 0:
            sys.exit("The TSV file is empty. Please check that the input VCF contains missense, inframe indel, or frameshift mutations.")
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

        if self.net_chop_method:
            self.net_chop()
        else:
            os.symlink(self.coverage_filter_out_path(), self.net_chop_out_path())
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
            + "Done: pvacseq has completed. File %s contains list of filtered putative neoantigens. " % self.final_path()
            + "We recommend appending coverage information and running `pvacseq coverage_filter` to filter based on sequencing coverage information"
        )
        if self.keep_tmp_files is False:
            shutil.rmtree(self.tmp_dir)

class MHCIPipeline(Pipeline):
    def __init__(self, **kwargs):
        Pipeline.__init__(self, **kwargs)
        self.peptide_sequence_length = kwargs['peptide_sequence_length']
        self.epitope_lengths         = kwargs['epitope_lengths']

    def generate_fasta(self, chunks):
        status_message("Generating Variant Peptide FASTA and Key Files")
        for (split_start, split_end) in chunks:
            tsv_chunk = "%d-%d" % (split_start, split_end)
            fasta_chunk = "%d-%d" % (split_start*2-1, split_end*2)
            split_tsv_file_path       = "%s_%s" % (self.tsv_file_path(), tsv_chunk)
            split_fasta_file_path     = "%s_%s" % (self.split_fasta_basename(), fasta_chunk)
            if os.path.exists(split_fasta_file_path):
                status_message("Split FASTA file for Entries %s already exists. Skipping." % (fasta_chunk))
                continue
            split_fasta_key_file_path = split_fasta_file_path + '.key'
            status_message("Generating Variant Peptide FASTA and Key Files - Entries %s" % (fasta_chunk))
            generate_fasta_params = [
                split_tsv_file_path,
                str(self.peptide_sequence_length),
                str(min(self.epitope_lengths)),
                split_fasta_file_path,
                split_fasta_key_file_path,
            ]
            if self.downstream_sequence_length:
                generate_fasta_params.extend(['-d', self.downstream_sequence_length,])
            lib.generate_fasta.main(generate_fasta_params)
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
                        lib.call_iedb.main([
                            split_fasta_file_path,
                            split_iedb_out,
                            iedb_method,
                            a,
                            '-l', str(epl),
                            '-r', str(self.iedb_retries),
                            '-e', self.iedb_executable,
                        ])
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
                        params = [
                            *split_iedb_output_files,
                            split_tsv_file_path,
                            split_fasta_key_file_path,
                            split_parsed_file_path,
                            '-m', self.top_score_metric,
                        ]
                        if self.top_result_per_mutation == True:
                            params.append('-t')
                        lib.parse_output.main(params)
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
            split_tsv_file_path       = "%s_%s" % (self.tsv_file_path(), tsv_chunk)
            split_fasta_file_path     = "%s_%s" % (self.split_fasta_basename(), fasta_chunk)
            if os.path.exists(split_fasta_file_path):
                status_message("Split FASTA file for Entries %s already exists. Skipping." % (fasta_chunk))
                continue
            split_fasta_key_file_path = split_fasta_file_path + '.key'
            status_message("Generating Variant Peptide FASTA and Key Files - Entries %s" % (fasta_chunk))
            generate_fasta_params = [
                split_tsv_file_path,
                str(self.peptide_sequence_length),
                '9', #This is the default core epitope length for IEDB class ii predictions
                split_fasta_file_path,
                split_fasta_key_file_path,
            ]
            if self.downstream_sequence_length:
                generate_fasta_params.extend(['-d', self.downstream_sequence_length,])
            lib.generate_fasta.main(generate_fasta_params)
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
                    lib.call_iedb.main([
                        split_fasta_file_path,
                        split_iedb_out,
                        iedb_method,
                        a,
                        '-r', str(self.iedb_retries),
                        '-e', self.iedb_executable,
                    ])
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
                    params = [
                        *split_iedb_output_files,
                        split_tsv_file_path,
                        split_fasta_key_file_path,
                        split_parsed_file_path,
                        '-m', self.top_score_metric,
                    ]
                    if self.top_result_per_mutation == True:
                        params.append('-t')
                    lib.parse_output.main(params)
                    status_message("Completed")
                    split_parsed_output_files.append(split_parsed_file_path)

        return split_parsed_output_files
