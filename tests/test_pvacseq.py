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

from pvactools.tools.pvacseq import *
import pvactools.tools.pvacseq.main as pvacseq_main
from .test_utils import *

def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacseq'
    )

class PvacseqTests(unittest.TestCase):
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

    def test_pvacseq_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacseq',
            "main.py"
        ))
        self.assertTrue(compiled_pvac_path)

    def test_parser(self):
        parser = pvacseq_main.define_parser()
        self.assertEqual(type(parser), argparse.ArgumentParser)

    def test_pvacseq_commands(self):
        pvac_script_path = os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacseq',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "allele_specific_cutoffs",
            "binding_filter",
            "coverage_filter",
            "download_example_data",
            "generate_aggregated_report",
            "generate_protein_fasta",
            "install_vep_plugin",
            "run",
            "net_chop",
            "netmhc_stab",
            "calculate_reference_proteome_similarity",
            "top_score_filter",
            "transcript_support_level_filter",
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
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "run.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_allele_specific_cutoffs_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
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
            "pvacseq",
            "binding_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_binding_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        binding_filter.main([input_file, output_file.name])

    def test_coverage_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "binding_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_coverage_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        coverage_filter.main([input_file, output_file.name])

    def test_download_example_data_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
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
            "pvacseq",
            "generate_aggregated_report.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_generate_aggregated_report_runs(self):
        input_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        generate_aggregated_report.main([input_file, output_file.name])

    def test_generate_protein_fasta_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "generate_protein_fasta.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_generate_protein_fasta_runs(self):
        input_file = os.path.join(self.test_data_directory, 'input.vcf')
        output_file = tempfile.NamedTemporaryFile()
        generate_protein_fasta.main([input_file, "25", output_file.name])
        os.unlink("{}.manufacturability.tsv".format(output_file.name))

    def test_install_vep_pugin_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "install_vep_plugin.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_install_vep_pugin_runs(self):
        output_dir = tempfile.TemporaryDirectory()
        install_vep_plugin.main([output_dir.name])

    def test_top_score_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "top_score_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_top_score_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        top_score_filter.main([input_file, output_file.name])

    def test_transcript_support_level_filter_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "transcript_support_level_filter.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_transcript_support_level_filter_runs(self):
        input_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'Test.all_epitopes.tsv')
        output_file = tempfile.NamedTemporaryFile()
        transcript_support_level_filter.main([input_file, output_file.name])

    def test_valid_alleles_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacseq",
            "valid_alleles.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_valid_alleles_runs(self):
        valid_alleles.main(["-p", "SMM"])

    def test_pvacseq_pipeline(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request, unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "input.vcf"),
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
                '--tdna-vaf', '20',
                '-d', 'full',
                '--pass-only',
                '--run-reference-proteome-similarity',
            ])

            run.main([
                os.path.join(self.test_data_directory, "input.vcf"),
                'sample.name',
                'DRB1*11:01',
                'NNalign',
                output_dir.name,
                '-e2', '15',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
                '-d', 'full',
                '--run-reference-proteome-similarity',
            ])
            close_mock_fhs()

            for file_name in (
                'sample.name.all_epitopes.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.tsv',
                'sample.name.tsv_1-24',
                'sample.name.fasta',
                'sample.name.all_epitopes.aggregated.tsv',
                'sample.name.filtered.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

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
                'sample.name.tsv',
                'sample.name.tsv_1-24',
                'sample.name.fasta',
                'sample.name.all_epitopes.aggregated.tsv',
                'sample.name.filtered.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', file_name.replace('sample.name', 'Test'))
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.all_epitopes.tsv',
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
                    os.path.join(self.test_data_directory, "input.vcf"),
                    'sample.name',
                    'DRB1*11:01',
                    'NNalign',
                    output_dir.name,
                    '-e2', '15',
                    '--top-score-metric=lowest',
                    '--keep-tmp-files',
                    '--run-reference-proteome-similarity',
                ])
            self.assertEqual(
                str(cm.exception),
                "Restart inputs are different from past inputs: \n" +
                "Past input: downstream_sequence_length - None\n" +
                "Current input: downstream_sequence_length - 1000\nAborting."
            )

            output_dir.cleanup()

    @patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
        data,
        files,
        test_data_directory()
    )))
    def test_pvacseq_pipeline_additional_report_columns(self):
        output_dir = tempfile.TemporaryDirectory()
        params = [
            os.path.join(self.test_data_directory, "input.vcf"),
            'Test',
            'DRB1*11:01',
            'NNalign',
            output_dir.name,
            '-e2', '15',
            '-a', 'sample_name',
            '-d', 'full',
            '--top-score-metric', 'lowest'
        ]
        run.main(params)
        output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'Test.filtered.tsv')
        expected_file = os.path.join(self.test_data_directory, 'Test_with_additional_report_columns.final.tsv')
        self.assertTrue(cmp(output_file, expected_file, False))
        output_dir.cleanup()

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
            '-e1', '8',
            '-s', '1000',
            '-k',
            '-p', os.path.join(self.test_data_directory, 'phased.vcf.gz'),
        ]
        run.main(params)

        for file_name in ['Test.8.fa.split_1-822', 'Test.8.fa.split_1-822.key']:
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
            expected_file = os.path.join(self.test_data_directory, 'phased', 'MHC_Class_I', 'tmp', file_name)
            self.assertTrue(cmp(output_file, expected_file, False), "files don't match %s - %s" %(output_file, expected_file))
        for file_name in ['Test.all_epitopes.tsv', 'Test.fasta']:
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
            expected_file = os.path.join(self.test_data_directory, 'phased', 'MHC_Class_I', file_name)
            self.assertTrue(compare(output_file, expected_file))
        for file_name in ['Test.proximal_variants.tsv', 'Test.all_epitopes.aggregated.tsv', 'Test.filtered.tsv']:
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
            expected_file = os.path.join(self.test_data_directory, 'phased', 'MHC_Class_I', file_name)
            self.assertTrue(cmp(output_file, expected_file, False), "files don't match %s - %s" %(output_file, expected_file))
        output_dir.cleanup()

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
            '-e1', '9,10',
            '-e2', '15',
            '--top-score-metric=lowest',
            '--keep-tmp-files',
            '-d', 'full',
        ])

        for file_name in (
            'Test.all_epitopes.tsv',
            'Test.filtered.tsv',
        ):
            output_file   = os.path.join(output_dir.name, 'combined', file_name)
            expected_file = os.path.join(self.test_data_directory, 'combine_and_condense', 'combined', file_name)
            self.assertTrue(compare(output_file, expected_file))
        output_dir.cleanup()

    def test_mismatched_allele_species_raises_exception(self):
        with self.assertRaises(Exception) as context:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
            run.main([
                os.path.join(self.test_data_directory, "input.vcf"),
                'Test',
                'HLA-G*01:09,H2-IAb',
                'NetMHC',
                output_dir.name,
                '-e1', '9',
            ])
            output_dir.cleanup()
        self.assertTrue('Requested alleles are not from the same species.' in str(context.exception))
