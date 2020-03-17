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
import lib
import datetime
from tools.pvacfuse import *
from mock import patch
from .test_utils import *
import tools.pvacfuse.main as pvacfuse_main
import argparse

def test_data_directory():
    return os.path.join(
        pvac_directory(),
        'tests',
        'test_data',
        'pvacfuse'
    )

class PvacfuseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pVac_directory = pvac_directory()
        cls.test_data_directory = test_data_directory()

    def test_pvacfuse_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pVac_directory,
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
            self.pVac_directory,
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
            "top_score_filter",
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
            self.pVac_directory,
            "tools",
            "pvacfuse",
            "run.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_allele_specific_cutoffs_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            "tools",
            "pvacfuse",
            "allele_specific_cutoffs.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_allele_specific_cutoffs_runs(self):
        allele_specific_cutoffs.main([])

    def test_binding_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pVac_directory,
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
            self.pVac_directory,
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
            self.pVac_directory,
            "tools",
            "pvacfuse",
            "top_score_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_top_score_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        top_score_filter.main([input_file, output_file.name])

    def test_valid_alleles_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pVac_directory,
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
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "fusions_annotated.bedpe"),
                'sample.name',
                'HLA-A*29:02',
                'NetMHC',
                output_dir.name,
                '-e', '9',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '--run-reference-proteome-similarity',
            ])

            for file_name in (
                'sample.name.tsv',
                'sample.name.tsv_1-4',
                'sample.name.fasta',
                'sample.name.all_epitopes.tsv',
                'sample.name.filtered.tsv',
                'sample.name.filtered.condensed.ranked.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file),  "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name_21.fa.split_1-8',
                'sample.name_21.fa.split_1-8.key',
                'sample.name.ann.HLA-A*29:02.9.tsv_1-8',
                'sample.name.HLA-A*29:02.9.parsed.tsv_1-8',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            mock_request.assert_has_calls([
                generate_class_i_call('ann', 'HLA-A*29:02', 9, os.path.join(output_dir.name, "MHC_Class_I", "tmp", "sample.name_21.fa.split_1-8"))
            ])

            output_dir.cleanup()

    def test_pvacfuse_combine_and_condense_steps(self):
        with unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
            for subdir in ['MHC_Class_I', 'MHC_Class_II']:
                path = os.path.join(output_dir.name, subdir)
                os.mkdir(path)
                test_data_dir = os.path.join(self.test_data_directory, 'combine_and_condense', subdir)
                for item in os.listdir(test_data_dir):
                    os.symlink(os.path.join(test_data_dir, item), os.path.join(path, item))

            run.main([
                os.path.join(self.test_data_directory, "fusions_annotated.bedpe"),
                'Test',
                'HLA-A*29:02,DRB1*11:01',
                'NetMHC', 'NNalign',
                output_dir.name,
                '-e', '9',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '--run-reference-proteome-similarity',
            ])
            close_mock_fhs()

            for file_name in (
                'Test.all_epitopes.tsv',
                'Test.filtered.tsv',
                'Test.filtered.condensed.ranked.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'combined', file_name)
                expected_file = os.path.join(self.test_data_directory, 'combine_and_condense', 'combined', file_name)
                self.assertTrue(compare(output_file, expected_file))
