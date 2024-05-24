import unittest
import unittest.mock
import os
import re
import sys
import tempfile
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run
from filecmp import cmp
import yaml
import datetime
from mock import patch
import argparse
from urllib.request import urlopen
from shutil import copyfileobj
from tempfile import NamedTemporaryFile

from pvactools.tools.pvacsplice import *
import pvactools.tools.pvacsplice.main as pvacsplice_main
from tests.utils import *

def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacsplice'
    )

class PvacspliceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pvactools_directory = pvactools_directory()
        cls.test_data_directory = test_data_directory()
        cls.methods = {
            'ann': {
                'HLA-E*01:01': [9, 10]
            },
            'pickpocket': {
                'HLA-G*01:09': [9, 10],
                'HLA-E*01:01': [9, 10],
            },
        }
        cls.peptide_fasta = os.path.join(pvactools_directory(), "tests", "test_data", "Homo_sapiens.GRCh38.pep.all.fa.gz")

    def test_pvacsplice_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacsplice',
            "main.py"
        ))
        self.assertTrue(compiled_pvac_path)

    def test_parser(self):
        parser = pvacsplice_main.define_parser()
        self.assertEqual(type(parser), argparse.ArgumentParser)

    def test_pvacsplice_commands(self):
        pvac_script_path = os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacsplice',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "allele_specific_cutoffs",
            "binding_filter",
            "coverage_filter",
            "run",
            "top_score_filter",
            "transcript_support_level_filter",
            "download_example_data",
            "generate_aggregated_report",
            "generate_protein_fasta",
            "install_vep_plugin",
            "net_chop",
            "netmhc_stab",
            "calculate_reference_proteome_similarity",
            "valid_alleles",
            'identify_problematic_amino_acids',
            ]:
            result = subprocess_run([
                sys.executable,
                pvac_script_path,
                command,
                '-h'
            ], shell=False, stdout=PIPE)
            self.assertFalse(result.returncode, "Failed `pvacsplice {} -h`".format(command))
            self.assertRegex(result.stdout.decode(), usage_search)

    def test_run_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "run.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_allele_specific_cutoffs_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "allele_specific_cutoffs.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_allele_specific_cutoffs_runs(self):
        allele_specific_cutoffs.main([])

    def test_binding_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "binding_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_binding_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'results', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        binding_filter.main([input_file, output_file.name])

    def test_coverage_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "coverage_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_coverage_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'results', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        coverage_filter.main([input_file, output_file.name])

    def test_transcript_support_level_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "transcript_support_level_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_transcript_support_level_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'results', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        transcript_support_level_filter.main([input_file, output_file.name])

    def test_top_score_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "top_score_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_top_score_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'results', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        top_score_filter.main([input_file, output_file.name])

    def test_download_example_data_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "download_example_data.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_download_example_data_runs(self):
        output_dir = tempfile.TemporaryDirectory()
        download_example_data.main([output_dir.name])

    def test_generate_aggregated_report_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "generate_aggregated_report.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_generate_aggregated_report_runs(self):
        input_file = os.path.join(self.test_data_directory, 'results', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        generate_aggregated_report.main([input_file, output_file.name])

    def test_install_vep_plugin_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "install_vep_plugin.py",
        ))
        self.assertTrue(compiled_run_path)

    def test_install_vep_plugin_runs(self):
        output_dir = tempfile.TemporaryDirectory()
        install_vep_plugin.main([output_dir.name])

    def test_net_chop_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "net_chop.py",
        ))
        self.assertTrue(compiled_run_path)

    def test_netmhc_stab_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "netmhc_stab.py",
        ))
        self.assertTrue(compiled_run_path)

    def test_calculate_reference_proteome_similarity_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "calculate_reference_proteome_similarity.py",
        ))
        self.assertTrue(compiled_run_path)

    def test_valid_alleles_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "valid_alleles.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_valid_alleles_runs(self):
        valid_alleles.main(["-p", "SMM"])

    def test_identify_problematic_amino_acids_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "identify_problematic_amino_acids.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_identify_problematic_amino_acids_runs(self):
        input_file = os.path.join(self.test_data_directory, 'results', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        identify_problematic_amino_acids.main([input_file, output_file.name, "C"])

    def test_pvacsplice_pipeline_class_I(self):
        with patch('pvactools.lib.call_iedb.requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            os.path.join(test_data_directory(), 'mock_files'),
        ))) as mock_request, patch('pvactools.lib.net_chop.NetChop.post_query', unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            os.path.join(test_data_directory(), 'mock_files'),
            'net_chop.html'
        ))), patch('pvactools.lib.netmhc_stab.NetMHCStab.query_netmhcstabpan_server',  unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            os.path.join(test_data_directory(), 'mock_files'),
            'Netmhcstab.html'
        ))):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            fasta_file = os.path.join(self.test_data_directory, "inputs", "all_sequences_chr1.fa.gz")
            unzipped_fasta_file = gunzip_file(fasta_file, suffix=".fa")
            run.main([
                os.path.join(self.test_data_directory, "inputs", "splice_junctions_chr1.tsv"),
                'HCC1395_TUMOR_DNA',
                'HLA-G*01:09,HLA-E*01:01',
                'NetMHC',
                'PickPocket',
                output_dir.name,
                os.path.join(self.test_data_directory, "inputs", "annotated.expression_chr1.vcf.gz"),
                unzipped_fasta_file,
                os.path.join(self.test_data_directory, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.gtf.gz"),
                '-e1', '9,10',
                '--normal-sample-name', 'HCC1395_NORMAL_DNA',
                '--keep-tmp-files',
                '-g',
                '--run-reference-proteome-similarity',
                '--peptide-fasta', self.peptide_fasta,
                '--net-chop-method', 'cterm',
                '--netmhc-stab',
                '--problematic-amino-acids', 'C',
                '-b', '2000',
                '--maximum-transcript-support-level', '3',
            ])

            close_mock_fhs()

            for file_name in (
                'HCC1395_TUMOR_DNA.transcripts.fa',
                'HCC1395_TUMOR_DNA_combined.tsv',
                'HCC1395_TUMOR_DNA_gtf.tsv',
            ):
                output_file   = os.path.join(output_dir.name, file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'HCC1395_TUMOR_DNA.all_epitopes.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run', 'MHC_Class_I', file_name)
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'HCC1395_TUMOR_DNA.all_epitopes.aggregated.tsv',
                'HCC1395_TUMOR_DNA.all_epitopes.aggregated.tsv.reference_matches',
                'HCC1395_TUMOR_DNA.filtered.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run', 'MHC_Class_I', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for length in [9, 10]:
                for file_name in (
                    'HCC1395_TUMOR_DNA.all_epitopes.tsv',
                ):
                    output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'MHC_Class_I_{}'.format(length), file_name)
                    expected_file = os.path.join(self.test_data_directory, 'results', 'run', 'MHC_Class_I', 'MHC_Class_I_{}'.format(length), file_name)
                    self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))
                for file_name in (
                    'HCC1395_TUMOR_DNA.{}.fa'.format(length),
                ):
                    output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'MHC_Class_I_{}'.format(length), 'tmp', file_name)
                    expected_file = os.path.join(self.test_data_directory, 'results', 'run', 'MHC_Class_I', 'MHC_Class_I_{}'.format(length), 'tmp', file_name)
                    self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            #Class I output files
            methods = self.methods
            for method in methods.keys():
                for allele in methods[method].keys():
                    mock_request.assert_has_calls([
                        generate_class_i_call(method, allele, 9, os.path.join(output_dir.name, "MHC_Class_I", "MHC_Class_I_9", "tmp", "HCC1395_TUMOR_DNA.9.fa.split_1-39"))
                    ])
                    mock_request.assert_has_calls([
                        generate_class_i_call(method, allele, 10, os.path.join(output_dir.name, "MHC_Class_I", "MHC_Class_I_10", "tmp", "HCC1395_TUMOR_DNA.10.fa.split_1-41"))
                    ])

            with self.assertRaises(SystemExit) as cm:
                run.main([
                    os.path.join(self.test_data_directory, "inputs", "splice_junctions_chr1.tsv"),
                    'HCC1395_TUMOR_DNA',
                    'HLA-G*01:09,HLA-E*01:01',
                    'NetMHC',
                    'PickPocket',
                    output_dir.name,
                    os.path.join(self.test_data_directory, "inputs", "annotated.expression_chr1.vcf.gz"),
                    unzipped_fasta_file,
                    os.path.join(self.test_data_directory, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.gtf.gz"),
                    '-e1', '9,10',
                    '--normal-sample-name', 'HCC1395_NORMAL_DNA',
                    '-b', '2000',
                    '--keep-tmp-files',
                    '-g',
                    '--run-reference-proteome-similarity',
                    '--peptide-fasta', self.peptide_fasta,
                    '--net-chop-method', 'cterm',
                    '--netmhc-stab',
                    '--problematic-amino-acids', 'C',
                    '-b', '2000',
                ])
            self.assertEqual(
                str(cm.exception),
                "Restart inputs are different from past inputs: \n" +
                "Past input: maximum_transcript_support_level - 3\n" +
                "Current input: maximum_transcript_support_level - 1\nAborting."
            )

            output_dir.cleanup()
            os.unlink(unzipped_fasta_file)

    def test_pvacsplice_pipeline_class_II(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            os.path.join(test_data_directory(), 'mock_files'),
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            fasta_file = os.path.join(self.test_data_directory, "inputs", "all_sequences_chr1.fa.gz")
            unzipped_fasta_file = gunzip_file(fasta_file, suffix=".fa")
            run.main([
                os.path.join(self.test_data_directory, "inputs", "splice_junctions_chr1.tsv"),
                'HCC1395_TUMOR_DNA',
                'DRB1*11:01',
                'NNalign',
                output_dir.name,
                os.path.join(self.test_data_directory, "inputs", "annotated.expression_chr1.vcf.gz"),
                unzipped_fasta_file,
                os.path.join(self.test_data_directory, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.gtf.gz"),
                '-e2', '15',
                '--normal-sample-name', 'HCC1395_NORMAL_DNA',
                '--keep-tmp-files',
                '-g',
                '--run-reference-proteome-similarity',
                '--peptide-fasta', self.peptide_fasta,
                '--problematic-amino-acids', 'C',
                '--maximum-transcript-support-level', '3',
                '-b', '2000',
            ])

            close_mock_fhs()

            #Class II output files
            for file_name in (
                'HCC1395_TUMOR_DNA.transcripts.fa',
                'HCC1395_TUMOR_DNA_combined.tsv',
                'HCC1395_TUMOR_DNA_gtf.tsv',
            ):
                output_file   = os.path.join(output_dir.name, file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'HCC1395_TUMOR_DNA.all_epitopes.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run', 'MHC_Class_II', file_name)
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'HCC1395_TUMOR_DNA.all_epitopes.aggregated.tsv',
                'HCC1395_TUMOR_DNA.all_epitopes.aggregated.tsv.reference_matches',
                'HCC1395_TUMOR_DNA.filtered.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run', 'MHC_Class_II', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'HCC1395_TUMOR_DNA.all_epitopes.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'MHC_Class_II_15', file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run', 'MHC_Class_II', 'MHC_Class_II_15', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'HCC1395_TUMOR_DNA.15.fa',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'MHC_Class_II_15', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run', 'MHC_Class_II', 'MHC_Class_II_15', 'tmp', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            mock_request.assert_has_calls([
                generate_class_ii_call('nn_align', 'DRB1*11:01', 15, os.path.join(output_dir.name, "MHC_Class_II", "MHC_Class_II_15", "tmp", "HCC1395_TUMOR_DNA.15.fa.split_1-51"))
            ])

            output_dir.cleanup()
            os.unlink(unzipped_fasta_file)

    def test_pvacsplice_combine_and_condense_steps(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            os.path.join(test_data_directory(), 'mock_files'),
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            fasta_file = os.path.join(self.test_data_directory, "inputs", "all_sequences_chr1.fa.gz")
            unzipped_fasta_file = gunzip_file(fasta_file, suffix=".fa")
            run.main([
                os.path.join(self.test_data_directory, "inputs", "splice_junctions_chr1.tsv"),
                'HCC1395_TUMOR_DNA',
                'HLA-G*01:09,HLA-E*01:01,DRB1*11:01',
                'NetMHC',
                'PickPocket',
                'NNalign',
                output_dir.name,
                os.path.join(self.test_data_directory, "inputs", "annotated.expression_chr1.vcf.gz"),
                unzipped_fasta_file,
                os.path.join(self.test_data_directory, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.gtf.gz"),
                '-e1', '9,10',
                '-e2', '15',
                '--normal-sample-name', 'HCC1395_NORMAL_DNA',
                '--keep-tmp-files',
                '-g',
                '--maximum-transcript-support-level', '3',
                '-b', '2000',
            ])

            for file_name in (
                'HCC1395_TUMOR_DNA.all_epitopes.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'combined', file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run_combined', file_name)
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'HCC1395_TUMOR_DNA.all_epitopes.aggregated.tsv',
                'HCC1395_TUMOR_DNA.filtered.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'combined', file_name)
                expected_file = os.path.join(self.test_data_directory, 'results', 'run_combined', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))
            output_dir.cleanup()
            os.unlink(unzipped_fasta_file)

    def test_mismatched_allele_species_raises_exception(self):
        with self.assertRaises(Exception) as context:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
            fasta_file = os.path.join(self.test_data_directory, "inputs", "all_sequences_chr1.fa.gz")
            unzipped_fasta_file = gunzip_file(fasta_file, suffix=".fa")
            run.main([
                os.path.join(self.test_data_directory, "inputs", "splice_junctions_chr1.tsv"),
                'HCC1395_TUMOR_DNA',
                'HLA-G*01:09,H2-IAb',
                'NetMHC',
                output_dir.name,
                os.path.join(self.test_data_directory, "inputs", "annotated.expression_chr1.vcf.gz"),
                unzipped_fasta_file,
                os.path.join(self.test_data_directory, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.gtf.gz"),
                '-e1', '9,10',
            ])
            output_dir.cleanup()
            os.unlink(unzipped_fasta_file)
        self.assertTrue('Requested alleles are not from the same species.' in str(context.exception))
