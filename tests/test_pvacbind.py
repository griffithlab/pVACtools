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

from pvactools.lib.pipeline import PvacbindPipeline
from pvactools.tools.pvacbind import *
from .test_utils import *

def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacbind'
    )

class PvacbindTests(unittest.TestCase):
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

    def test_pvacbind_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            'tools',
            'pvacbind',
            "main.py"
        ))
        self.assertTrue(compiled_pvac_path)

    def test_pvacbind_commands(self):
        pvac_script_path = os.path.join(
            self.pvactools_directory,
            "pvactools",
            'tools',
            'pvacbind',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "run",
            'binding_filter',
            'valid_alleles',
            'allele_specific_cutoffs',
            'download_example_data',
            "net_chop",
            "netmhc_stab",
            "calculate_reference_proteome_similarity",
            'top_score_filter',
            'generate_aggregated_report',
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
            "pvacbind",
            "run.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_generate_aggregated_report_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            "tools",
            "pvacbind",
            "generate_aggregated_report.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_generate_aggregated_report_runs(self):
        input_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        generate_aggregated_report.main([input_file, output_file.name])

    def test_process_stops(self):
        output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
        params = {
            'input_file': os.path.join(self.test_data_directory, "input_with_stops.fasta"),
            'input_file_type': 'fasta',
            'sample_name': 'Test',
            'alleles': ['HLA-G*01:09'],
            'prediction_algorithms': ['NetMHC'],
            'output_dir': output_dir.name,
            'epitope_lengths': [9],
        }
        pipeline = PvacbindPipeline(**params)
        pipeline.create_per_length_fasta_and_process_stops(9)
        output_file   = os.path.join(output_dir.name, 'tmp', 'Test.9.fa')
        expected_file = os.path.join(self.test_data_directory, 'output_with_stops.fasta')
        self.assertTrue(cmp(output_file, expected_file))
        output_dir.cleanup()

    def test_pvacbind_pipeline(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, timeout=None, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request, unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "input.fasta"),
                'sample.name',
                'HLA-G*01:09,HLA-E*01:01',
                'NetMHC',
                'PickPocket',
                output_dir.name,
                '-e1', '9,10',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '--net-chop-method', 'cterm',
                '--netmhc-stab',
                '--run-reference-proteome-similarity',
            ])

            run.main([
                os.path.join(self.test_data_directory, "input.fasta"),
                'sample.name',
                'DRB1*11:01',
                'NNalign',
                output_dir.name,
                '-e2', '15',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '--run-reference-proteome-similarity',
            ])
            close_mock_fhs()

            for file_name in (
                'sample.name.all_epitopes.tsv',
                'sample.name.filtered.tsv',
                'sample.name.all_epitopes.aggregated.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.9.fa.split_1-48',
                'sample.name.9.fa.split_1-48.key',
                'sample.name.10.fa.split_1-48',
                'sample.name.10.fa.split_1-48.key',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.HLA-G*01:09.9.parsed.tsv_1-48',
                'sample.name.HLA-G*01:09.10.parsed.tsv_1-48',
                'sample.name.HLA-E*01:01.9.parsed.tsv_1-48',
                'sample.name.HLA-E*01:01.10.parsed.tsv_1-48',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

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
                            generate_class_i_call(method, allele, length, os.path.join(output_dir.name, "MHC_Class_I", "tmp", "sample.name.{}.fa.split_1-48".format(length)))
                        ])
                        output_file   = os.path.join(output_dir.name, "MHC_Class_I", "tmp", 'sample.name.%s.%s.%s.tsv_1-48' % (method, allele, length))
                        expected_file = os.path.join(self.test_data_directory, "MHC_Class_I", "tmp", 'Test.%s.%s.%s.tsv_1-48' % (method, allele, length))
                        self.assertTrue(cmp(output_file, expected_file, False), "files don't match %s - %s" %(output_file, expected_file))

            #Class II output files
            for file_name in (
                'sample.name.all_epitopes.tsv',
                'sample.name.filtered.tsv',
                'sample.name.all_epitopes.aggregated.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.15.fa.split_1-48',
                'sample.name.15.fa.split_1-48.key',
                'sample.name.nn_align.DRB1*11:01.15.tsv_1-48',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(cmp(output_file, expected_file, False), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.DRB1*11:01.15.parsed.tsv_1-48',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'inputs.yml',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'log', file_name)
                self.assertTrue(os.path.exists(output_file))

            mock_request.assert_has_calls([
                generate_class_ii_call('nn_align', 'DRB1*11:01', 15, os.path.join(output_dir.name, "MHC_Class_II", "tmp", "sample.name.15.fa.split_1-48"))
            ])

            with self.assertRaises(SystemExit) as cm:
                run.main([
                    os.path.join(self.test_data_directory, "input.fasta"),
                    'sample.name',
                    'DRB1*11:01',
                    'NNalign',
                    output_dir.name,
                    '-e2', '15',
                    '--keep-tmp-files',
                    '--run-reference-proteome-similarity',
                ])
            self.assertEqual(
                str(cm.exception),
                "Restart inputs are different from past inputs: \n" +
                "Past input: top_score_metric - lowest\n" +
                "Current input: top_score_metric - median\nAborting."
            )

            output_dir.cleanup()

    def test_duplicate_fasta_header(self):
        with self.assertRaises(Exception) as cm:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
            run.main([
                os.path.join(self.test_data_directory, "input.duplicate_header.fasta"),
                'Test',
                'HLA-A*02:01',
                'NetMHC',
                output_dir.name,
                '-e1', '8'
            ])
        self.assertEqual(
            str(cm.exception),
            "Duplicate fasta header 1. Please ensure that the input FASTA uses unique headers."
        )
        output_dir.cleanup()

    def test_pvacbind_combine_and_condense_steps(self):
        with unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
            for subdir in ['MHC_Class_I', 'MHC_Class_II']:
                path = os.path.join(output_dir.name, subdir)
                os.mkdir(path)
                test_data_dir = os.path.join(self.test_data_directory, 'combine_and_condense', subdir)
                for item in os.listdir(test_data_dir):
                    os.symlink(os.path.join(test_data_dir, item), os.path.join(path, item))

            run.main([
                os.path.join(self.test_data_directory, "input.fasta"),
                'Test',
                'HLA-G*01:09,HLA-E*01:01,DRB1*11:01',
                'NetMHC',
                'PickPocket',
                'NNalign',
                output_dir.name,
                '-e1', '9,10',
                '-e2', '15',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '--allele-specific-binding-thresholds',
                '--run-reference-proteome-similarity',
            ])
            close_mock_fhs()

            for file_name in (
                'Test.all_epitopes.tsv',
                'Test.filtered.tsv',
                'Test.all_epitopes.aggregated.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'combined', file_name)
                expected_file = os.path.join(self.test_data_directory, 'combine_and_condense', 'combined', file_name)
                self.assertTrue(compare(output_file, expected_file))
