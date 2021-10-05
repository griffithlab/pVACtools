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

from pvactools.tools.pvacfuse import *
import pvactools.tools.pvacfuse.main as pvacfuse_main
from .test_utils import *

def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacfuse'
    )

class PvacfuseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pvactools_directory = pvactools_directory()
        cls.test_data_directory = test_data_directory()

    def test_pvacfuse_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            'tools',
            'pvacfuse',
            'main.py'
        ))
        self.assertTrue(compiled_pvac_path)

    def test_parser(self):
        parser = pvacfuse_main.define_parser()
        self.assertEqual(type(parser), argparse.ArgumentParser)

    def test_pvacfuse_commands(self):
        pvac_script_path = os.path.join(
            self.pvactools_directory,
            "pvactools",
            'tools',
            'pvacfuse',
            'main.py'
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "run",
            "allele_specific_cutoffs",
            "binding_filter",
            "valid_alleles",
            "download_example_data",
            "net_chop",
            "netmhc_stab",
            "calculate_reference_proteome_similarity",
            "top_score_filter",
            "generate_aggregated_report",
            ]:
            result = subprocess_run([
                sys.executable,
                pvac_script_path,
                command,
                '-h'
            ], shell=False, stdout=PIPE)
            self.assertFalse(result.returncode)
            self.assertRegex(result.stdout.decode(), usage_search)

    def test_run_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            "tools",
            "pvacfuse",
            "run.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_allele_specific_cutoffs_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            "tools",
            "pvacfuse",
            "allele_specific_cutoffs.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_allele_specific_cutoffs_runs(self):
        allele_specific_cutoffs.main([])

    def test_binding_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            "tools",
            "pvacfuse",
            "binding_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_binding_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        binding_filter.main([input_file, output_file.name])

    def test_download_example_data_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            "tools",
            "pvacfuse",
            "download_example_data.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_download_example_data_runs(self):
        output_dir = tempfile.TemporaryDirectory()
        download_example_data.main([output_dir.name])

    def test_top_score_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            "tools",
            "pvacfuse",
            "top_score_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_top_score_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        top_score_filter.main([input_file, output_file.name])

    def test_generate_aggregated_report_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            "tools",
            "pvacfuse",
            "generate_aggregated_report.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_generate_aggregated_report_runs(self):
        input_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        generate_aggregated_report.main([input_file, output_file.name])

    def test_valid_alleles_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            "tools",
            "pvacfuse",
            "valid_alleles.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_valid_alleles_runs(self):
        valid_alleles.main(["-p", "SMM"])

    def test_pvacfuse_pipeline(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request, unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "agfusion"),
                'sample.name',
                'HLA-A*29:02',
                'NetMHC',
                output_dir.name,
                '-e1', '9',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '--run-reference-proteome-similarity',
            ])
            close_mock_fhs()

            for file_name in (
                'sample.name.fasta',
                'sample.name.all_epitopes.tsv',
                'sample.name.filtered.tsv',
                'sample.name.all_epitopes.aggregated.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file),  "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.ann.HLA-A*29:02.9.tsv_1-44',
                'sample.name.HLA-A*29:02.9.parsed.tsv_1-44',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', '9', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', '9', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            mock_request.assert_has_calls([
                generate_class_i_call('ann', 'HLA-A*29:02', 9, os.path.join(output_dir.name, "MHC_Class_I", "9", "tmp", "sample.name.9.fa.split_1-44"))
            ])

            output_dir.cleanup()

    def test_pvacfuse_combine_and_condense_steps(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request, unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
            run.main([
                os.path.join(self.test_data_directory, "agfusion"),
                'Test',
                'HLA-A*29:02,DRB1*11:01',
                'NetMHC', 'NNalign',
                output_dir.name,
                '-e1', '9',
                '-e2', '15',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '--run-reference-proteome-similarity',
            ])
            close_mock_fhs()

            for file_name in (
                'Test.all_epitopes.tsv',
                'Test.filtered.tsv',
                'Test.all_epitopes.aggregated.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'combined', file_name)
                expected_file = os.path.join(self.test_data_directory, 'combined', file_name)
                self.assertTrue(compare(output_file, expected_file))
