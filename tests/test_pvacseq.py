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
from lib.pipeline import *
import datetime
from tools.pvacseq import *
from mock import patch
from .test_utils import *

def test_data_directory():
    return os.path.join(
        pvac_directory(),
        'tests',
        'test_data',
        'pvacseq'
    )

class PvacseqTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pVac_directory = pvac_directory()
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

    def test_pvacseq_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            'tools',
            'pvacseq',
            "main.py"
        ))
        self.assertTrue(compiled_pvac_path)

    def test_pvacseq_commands(self):
        pvac_script_path = os.path.join(
            self.pVac_directory,
            'tools',
            'pvacseq',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "binding_filter",
            "coverage_filter",
            "run",
            "generate_protein_fasta",
            "install_vep_plugin",
            "download_example_data",
            "valid_alleles",
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
            "pvacseq",
            "run.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_pvacseq_pipeline(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request, unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "input.vcf"),
                'Test',
                'HLA-G*01:09,HLA-E*01:01',
                'NetMHC',
                'PickPocket',
                output_dir.name,
                '-e', '9,10',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '--net-chop-method', 'cterm',
                '--netmhc-stab',
                '--tdna-vaf', '20',
                '-d', 'full',
                '--pass-only',
            ])

            run.main([
                os.path.join(self.test_data_directory, "input.vcf"),
                'Test',
                'DRB1*11:01',
                'NNalign',
                output_dir.name,
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '-d', 'full',
            ])
            close_mock_fhs()

            for file_name in (
                'Test.all_epitopes.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', file_name)
                self.assertTrue(compare(output_file, expected_file))

            for file_name in (
                'Test.tsv',
                'Test.tsv_1-24',
                'Test.fasta',
                'Test.filtered.tsv',
                'Test.filtered.condensed.ranked.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', file_name)
                self.assertTrue(cmp(output_file, expected_file))

            for file_name in (
                'Test_21.fa.split_1-48',
                'Test_21.fa.split_1-48.key',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'tmp', file_name)
                self.assertTrue(cmp(output_file, expected_file))

            for file_name in (
                'Test.HLA-G*01:09.9.parsed.tsv_1-48',
                'Test.HLA-G*01:09.10.parsed.tsv_1-48',
                'Test.HLA-E*01:01.9.parsed.tsv_1-48',
                'Test.HLA-E*01:01.10.parsed.tsv_1-48',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'tmp', file_name)
                self.assertTrue(compare(output_file, expected_file))

            for file_name in (
                'inputs.yml',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'log', file_name)
                self.assertTrue(os.path.exists(output_file))

            #Class I output files
            methods = self.methods
            for method in methods.keys():
                for allele in methods[method].keys():
                    for length in methods[method][allele]:
                        mock_request.assert_has_calls([
                            generate_class_i_call(method, allele, length, os.path.join(output_dir.name, "MHC_Class_I", "tmp", "Test_21.fa.split_1-48"))
                        ])
                        output_file   = os.path.join(output_dir.name, "MHC_Class_I", "tmp", 'Test.%s.%s.%s.tsv_1-48' % (method, allele, length))
                        expected_file = os.path.join(self.test_data_directory, "MHC_Class_I", "tmp", 'Test.%s.%s.%s.tsv_1-48' % (method, allele, length))
                        self.assertTrue(cmp(output_file, expected_file, False))

            #Class II output files
            for file_name in (
                'Test.tsv',
                'Test.tsv_1-24',
                'Test.fasta',
                'Test.all_epitopes.tsv',
                'Test.filtered.tsv',
                'Test.filtered.condensed.ranked.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', file_name)
                self.assertTrue(compare(output_file, expected_file))

            for file_name in (
                'Test_31.fa.split_1-48',
                'Test_31.fa.split_1-48.key',
                'Test.nn_align.DRB1*11:01.15.tsv_1-48',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', 'tmp', file_name)
                self.assertTrue(cmp(output_file, expected_file, False))

            for file_name in (
                'Test.DRB1*11:01.15.parsed.tsv_1-48',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', 'tmp', file_name)
                self.assertTrue(compare(output_file, expected_file))

            for file_name in (
                'inputs.yml',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'log', file_name)
                self.assertTrue(os.path.exists(output_file))

            mock_request.assert_has_calls([
                generate_class_ii_call('nn_align', 'DRB1*11:01', os.path.join(output_dir.name, "MHC_Class_II", "tmp", "Test_31.fa.split_1-48"))
            ])

            with self.assertRaises(SystemExit) as cm:
                run.main([
                    os.path.join(self.test_data_directory, "input.vcf"),
                    'Test',
                    'DRB1*11:01',
                    'NNalign',
                    output_dir.name,
                    '--top-score-metric=lowest',
                    '--keep-tmp-files',
                ])
                self.assertEqual(
                    cm.exception,
                    "Restart inputs are different from past inputs: \n" +
                    "Past input: downstream_sequence_length - None\n" +
                    "Current input: downstream_sequence_length - 1000"
                )

            output_dir.cleanup()

    @patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
        data,
        files,
        test_data_directory()
    )))
    def test_pvacseq_pipeline_additional_report_columns(self):
        with unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            output_dir = tempfile.TemporaryDirectory()
            params = [
                os.path.join(self.test_data_directory, "input.vcf"),
                'Test',
                'HLA-E*01:01',
                'NetMHC',
                output_dir.name,
                '-e', '9,10',
                '-a', 'sample_name',
                '-d', 'full'
            ]
            run.main(params)
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'Test.filtered.tsv')
            expected_file = os.path.join(self.test_data_directory, 'Test_with_additional_report_columns.final.tsv')
            self.assertTrue(cmp(output_file, expected_file, False))
            close_mock_fhs()

    @patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
        data,
        files,
        test_data_directory()
    )))
    def test_pvacseq_pipeline_proximal_variants_vcf(self):
        output_dir = tempfile.TemporaryDirectory()

        params = [
            os.path.join(self.test_data_directory, "input_somatic.vcf.gz"),
            'Test',
            'HLA-E*01:01',
            'NetMHC',
            output_dir.name,
            '-e', '8',
            '-s', '1000',
            '-k',
            '-p', os.path.join(self.test_data_directory, 'phased.vcf.gz')
        ]
        run.main(params)

        for file_name in ['Test_21.fa.split_1-818', 'Test_21.fa.split_1-818.key']:
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
            expected_file = os.path.join(self.test_data_directory, 'phased', 'MHC_Class_I', 'tmp', file_name)
            self.assertTrue(cmp(output_file, expected_file, False))
        for file_name in ['Test.all_epitopes.tsv']:
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
            expected_file = os.path.join(self.test_data_directory, 'phased', 'MHC_Class_I', file_name)
            self.assertTrue(compare(output_file, expected_file))
        for file_name in ['Test.proximal_variants.tsv', 'Test.filtered.tsv']:
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
            expected_file = os.path.join(self.test_data_directory, 'phased', 'MHC_Class_I', file_name)
            self.assertTrue(cmp(output_file, expected_file, False))

    def test_pvacseq_combine_and_condense_steps(self):
        output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
        for subdir in ['MHC_Class_I', 'MHC_Class_II']:
            path = os.path.join(output_dir.name, subdir)
            os.mkdir(path)
            test_data_dir = os.path.join(self.test_data_directory, 'combine_and_condense', subdir)
            for item in os.listdir(test_data_dir):
                os.symlink(os.path.join(test_data_dir, item), os.path.join(path, item))

        run.main([
            os.path.join(self.test_data_directory, "input.vcf"),
            'Test',
            'HLA-G*01:09,HLA-E*01:01,DRB1*11:01',
            'NetMHC',
            'PickPocket',
            'NNalign',
            output_dir.name,
            '-e', '9,10',
            '--top-score-metric=lowest',
            '--keep-tmp-files',
            '--tdna-vaf', '20',
            '-d', 'full',
        ])

        for file_name in (
            'Test.all_epitopes.tsv',
            'Test.filtered.tsv',
            'Test.filtered.condensed.ranked.tsv',
        ):
            output_file   = os.path.join(output_dir.name, 'combined', file_name)
            expected_file = os.path.join(self.test_data_directory, 'combine_and_condense', 'combined', file_name)
            self.assertTrue(compare(output_file, expected_file))
