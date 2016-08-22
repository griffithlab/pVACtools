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

class Pipeline(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        self.input_file                  = kwargs['input_file']
        self.sample_name                 = kwargs['sample_name']
        self.alleles                     = kwargs['alleles']
        self.prediction_algorithms       = kwargs['prediction_algorithms']
        self.output_dir                  = kwargs['output_dir']
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
        self.downstream_sequence_length  = kwargs['downstream_sequence_length']
        self.keep_tmp_files              = kwargs['keep_tmp_files']
        tmp_dir = os.path.join(self.output_dir, 'tmp')
        os.makedirs(tmp_dir, exist_ok=True)
        self.tmp_dir = tmp_dir

    def tsv_file_path(self):
        tsv_file = self.sample_name + '.tsv'
        return os.path.join(self.output_dir, tsv_file)

    def convert_vcf(self):
        print("Converting VCF to TSV")
        sys.stdout.flush()
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
        print("Completed")

    def fasta_file_path(self):
        fasta_file = self.sample_name + "_" + str(self.peptide_sequence_length) + ".fa"
        return os.path.join(self.output_dir, fasta_file)

    @abstractmethod
    def generate_fasta(self):
        pass

    def fasta_entry_count(self):
        with open(self.fasta_file_path()) as f:
            for i, l in enumerate(f, 1):
                pass
        return i/2

    def split_fasta_basename(self):
        return os.path.join(self.tmp_dir, self.sample_name + "_" + str(self.peptide_sequence_length) + ".fa.split")

    def split_fasta_file_and_create_key_files(self):
        entry_count = self.fasta_entry_count()
        split_reader = open(self.fasta_file_path(), mode='r')
        split_start = 1
        #Each fasta entry consists of two lines: header and sequence
        chunk_size  = self.fasta_size * 2
        chunks = []
        for chunk in split_file(split_reader, chunk_size):
            split_end = split_start + self.fasta_size - 1
            if split_end > entry_count:
                split_end = entry_count
            print("Splitting FASTA into smaller chunks - Entries %d-%d" % (split_start, split_end))

            split_fasta_file_path = "%s_%d-%d"%(self.split_fasta_basename(), split_start, split_end)
            if os.path.exists(split_fasta_file_path):
                print("Split FASTA file for Entries %d-%d already exists. Skipping." % (split_start, split_end))

                [entry for entry in chunk]
            else:
                split_writer = open(split_fasta_file_path, mode='w')
                split_writer.writelines(chunk)
                split_writer.close()
                print("Completed")
            print("Generating FASTA Key File - Entries %d-%d" % (split_start, split_end))

            split_fasta_key_file_path = split_fasta_file_path + '.key'
            if os.path.exists(split_fasta_key_file_path):
                print("Split FASTA Key File for Entries %d-%d already exists. Skipping." % (split_start, split_end))

            else:
                lib.generate_fasta_key.main([
                    split_fasta_file_path,
                    split_fasta_key_file_path,
                ])
                print("Completed")
            chunks.append("%d-%d" % (split_start, split_end))
            split_start += self.fasta_size
        split_reader.close()
        return chunks

    @abstractmethod
    def call_iedb_and_parse_outputs(self, chunks):
        pass

    def combined_parsed_path(self):
        combined_parsed = "%s.combined.parsed.tsv" % self.sample_name
        return os.path.join(self.output_dir, combined_parsed)

    def combined_parsed_outputs(self, split_parsed_output_files):
        print("Combining Parsed IEDB Output Files")
        sys.stdout.flush()
        lib.combine_parsed_outputs.main([
            *split_parsed_output_files,
            self.combined_parsed_path()
        ])
        print("Completed")

    def binding_filter_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".filtered.binding.tsv")

    def binding_filter(self):
        print("Running Binding Filters")
        sys.stdout.flush()
        lib.binding_filter.main(
            [
                self.combined_parsed_path(),
                self.binding_filter_out_path(),
                '-c', str(self.minimum_fold_change),
                '-b', str(self.binding_threshold),
                '-m', str(self.top_score_metric),
            ]
        )
        print("Completed")

    def coverage_filter_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".filtered.coverage.tsv")

    def coverage_filter(self):
        print("Running Coverage Filters")
        sys.stdout.flush()
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
        print("Completed")

    def net_chop_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".chop.tsv")

    def net_chop(self):
        print("Submitting remaining epitopes to NetChop")
        sys.stdout.flush()
        lib.net_chop.main([
            self.coverage_filter_out_path(),
            self.net_chop_out_path(),
            '--method',
            self.net_chop_method,
            '--threshold',
            str(self.net_chop_threshold)
        ])
        print("Completed")

    def netmhc_stab_out_path(self):
        return os.path.join(self.output_dir, self.sample_name+".stab.tsv")

    def call_netmhc_stab(self):
        print("Running NetMHCStabPan")
        sys.stdout.flush()
        lib.netmhc_stab.main([
            self.net_chop_out_path(),
            self.netmhc_stab_out_path(),
        ])
        print("Completed")

    def final_path(self):
        return os.path.join(self.output_dir, self.sample_name+".final.tsv")

    def execute(self):
        self.convert_vcf()
        self.generate_fasta()

        if os.path.getsize(self.fasta_file_path()) == 0:
            sys.exit("The fasta file is empty. Please check that the input VCF contains missense, inframe indel, or frameshift mutations.")

        chunks                    = self.split_fasta_file_and_create_key_files()
        split_parsed_output_files = self.call_iedb_and_parse_outputs(chunks)

        if len(split_parsed_output_files) == 0:
            print("No output files were created. Aborting.")
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


        print("\n")
        print("Done: pvacseq has completed. File %s contains list of filtered putative neoantigens" % self.final_path())
        print("We recommend appending coverage information and running `pvacseq coverage_filter` to filter based on sequencing coverage information")
        sys.stdout.flush()
        if self.keep_tmp_files is False:
            shutil.rmtree(self.tmp_dir)

class MHCIPipeline(Pipeline):
    def __init__(self, **kwargs):
        Pipeline.__init__(self, **kwargs)
        self.peptide_sequence_length = kwargs['peptide_sequence_length']
        self.epitope_lengths         = kwargs['epitope_lengths']

    def generate_fasta(self):
        print("Generating Variant Peptide FASTA File")
        sys.stdout.flush()
        generate_fasta_params = [
            self.tsv_file_path(),
            str(self.peptide_sequence_length),
            str(min(self.epitope_lengths)),
            self.fasta_file_path(),
        ]
        if self.downstream_sequence_length:
            generate_fasta_params.extend(['-d', self.downstream_sequence_length,])
        lib.generate_fasta.main(generate_fasta_params)
        print("Completed")

    def call_iedb_and_parse_outputs(self, chunks):
        split_parsed_output_files = []
        for chunk in chunks:
            for a in self.alleles:
                for epl in self.epitope_lengths:
                    split_fasta_file_path = "%s_%s"%(self.split_fasta_basename(), chunk)
                    split_iedb_output_files = []
                    print("Processing entries for Allele %s and Epitope Length %s - Entries %s" % (a, epl, chunk))
                    for method in self.prediction_algorithms:
                        prediction_class = globals()[method]
                        prediction = prediction_class()
                        iedb_method = prediction.iedb_prediction_method
                        valid_alleles = prediction.valid_allele_names()
                        if a not in valid_alleles:
                            print("Allele %s not valid for Method %s. Skipping." % (a, method))
                            continue
                        valid_lengths = prediction.valid_lengths_for_allele(a)
                        if epl not in valid_lengths:
                            print("Epitope Length %s is not valid for Method %s and Allele %s. Skipping." % (epl, method, a))
                            continue

                        split_iedb_out = os.path.join(self.tmp_dir, ".".join([self.sample_name, iedb_method, a, str(epl), "tsv_%s" % chunk]))
                        if os.path.exists(split_iedb_out):
                            print("IEDB file for Allele %s and Epitope Length %s with Method %s (Entries %s) already exists. Skipping." % (a, epl, method, chunk))
                            split_iedb_output_files.append(split_iedb_out)
                            continue
                        print("Running IEDB on Allele %s and Epitope Length %s with Method %s - Entries %s" % (a, epl, method, chunk))
                        sys.stdout.flush()
                        lib.call_iedb.main([
                            split_fasta_file_path,
                            split_iedb_out,
                            iedb_method,
                            a,
                            '-l', str(epl),
                        ])
                        print("Completed")
                        split_iedb_output_files.append(split_iedb_out)

                    split_parsed_file_path = os.path.join(self.tmp_dir, ".".join([self.sample_name, a, str(epl), "parsed", "tsv_%s" % chunk]))
                    if os.path.exists(split_parsed_file_path):
                        print("Parsed Output File for Allele %s and Epitope Length %s (Entries %s) already exists. Skipping" % (a, epl, chunk))
                        split_parsed_output_files.append(split_parsed_file_path)
                        continue
                    split_fasta_key_file_path = split_fasta_file_path + '.key'

                    if len(split_iedb_output_files) > 0:
                        print("Parsing IEDB Output for Allele %s and Epitope Length %s - Entries %s" % (a, epl, chunk))
                        params = [
                            *split_iedb_output_files,
                            self.tsv_file_path(),
                            split_fasta_key_file_path,
                            split_parsed_file_path,
                            '-m', self.top_score_metric,
                        ]
                        if self.top_result_per_mutation == True:
                            params.append('-t')
                        lib.parse_output.main(params)
                        print("Completed")
                        split_parsed_output_files.append(split_parsed_file_path)
        return split_parsed_output_files

class MHCIIPipeline(Pipeline):
    def __init__(self, **kwargs):
        Pipeline.__init__(self, **kwargs)
        self.peptide_sequence_length = 31

    def generate_fasta(self):
        print("Generating Variant Peptide FASTA File")
        sys.stdout.flush()
        generate_fasta_params = [
            self.tsv_file_path(),
            str(self.peptide_sequence_length),
            '9', #This is the default core epitope length for IEDB class ii predictions
            self.fasta_file_path(),
        ]
        if self.downstream_sequence_length:
            generate_fasta_params.extend(['-d', self.downstream_sequence_length,])
        lib.generate_fasta.main(generate_fasta_params)
        print("Completed")

    def call_iedb_and_parse_outputs(self, chunks):
        split_parsed_output_files = []
        for chunk in chunks:
            for a in self.alleles:
                split_fasta_file_path = "%s_%s"%(self.split_fasta_basename(), chunk)
                split_iedb_output_files = []
                print("Processing entries for Allele %s - Entries %s" % (a, chunk))
                for method in self.prediction_algorithms:
                    prediction_class = globals()[method]
                    prediction = prediction_class()
                    iedb_method = prediction.iedb_prediction_method
                    valid_alleles = prediction.valid_allele_names()
                    if a not in valid_alleles:
                        print("Allele %s not valid for Method %s. Skipping." % (a, method))
                        continue

                    split_iedb_out = os.path.join(self.tmp_dir, ".".join([self.sample_name, iedb_method, a, "tsv_%s" % chunk]))
                    if os.path.exists(split_iedb_out):
                        print("IEDB file for Allele %s with Method %s (Entries %s) already exists. Skipping." % (a, method, chunk))
                        split_iedb_output_files.append(split_iedb_out)
                        continue
                    print("Running IEDB on Allele %s with Method %s - Entries %s" % (a, method, chunk))
                    sys.stdout.flush()
                    lib.call_iedb.main([
                        split_fasta_file_path,
                        split_iedb_out,
                        iedb_method,
                        a,
                    ])
                    print("Completed")
                    split_iedb_output_files.append(split_iedb_out)

                split_parsed_file_path = os.path.join(self.tmp_dir, ".".join([self.sample_name, a, "parsed", "tsv_%s" % chunk]))
                if os.path.exists(split_parsed_file_path):
                    print("Parsed Output File for Allele %s (Entries %s) already exists. Skipping" % (a, chunk))
                    split_parsed_output_files.append(split_parsed_file_path)
                    continue
                split_fasta_key_file_path = split_fasta_file_path + '.key'

                if len(split_iedb_output_files) > 0:
                    print("Parsing IEDB Output for Allele %s - Entries %s" % (a, chunk))
                    params = [
                        *split_iedb_output_files,
                        self.tsv_file_path(),
                        split_fasta_key_file_path,
                        split_parsed_file_path,
                        '-m', self.top_score_metric,
                    ]
                    if self.top_result_per_mutation == True:
                        params.append('-t')
                    lib.parse_output.main(params)
                    print("Completed")
                    split_parsed_output_files.append(split_parsed_file_path)

        return split_parsed_output_files


def split_file(reader, lines=400):
    from itertools import islice, chain
    tmp = next(reader)
    while tmp!="":
        yield chain([tmp], islice(reader, lines-1))
        try:
            tmp = next(reader)
        except StopIteration:
            return
