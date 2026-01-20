import unittest
import unittest.mock
from unittest.mock import patch
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
import argparse
from urllib.request import urlopen
from shutil import copyfileobj, copyfile
from tempfile import NamedTemporaryFile

from pvactools.tools.pvacseq import *
import pvactools.tools.pvacseq.main as pvacseq_main
from tests.utils import *

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
        cls.peptide_fasta = os.path.join(pvactools_directory(), "tests", "test_data", "Homo_sapiens.GRCh38.pep.short.fa.gz")

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
            "run",
            ]:
            result = subprocess_run([
                sys.executable,
                pvac_script_path,
                command,
                '-h'
            ], shell=False, stdout=PIPE)
            self.assertFalse(result.returncode, "Failed `pvacseq {} -h`".format(command))
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

    def test_pvacseq_pipeline(self):
        with patch('pvactools.lib.call_iedb.requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request, patch('pvactools.lib.net_chop.NetChop.post_query', unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            'net_chop.html'
        ))), patch('pvactools.lib.netmhc_stab.NetMHCStab.query_netmhcstabpan_server',  unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            'Netmhcstab.html'
        ))):
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
                '--top-score-metric2=ic50',
                '--keep-tmp-files',
                '--net-chop-method', 'cterm',
                '--netmhc-stab',
                '--tdna-vaf', '0.2',
                '-d', 'full',
                '--pass-only',
                '--run-reference-proteome-similarity',
                '--peptide-fasta', self.peptide_fasta,
                '--biotypes', 'IG_V_gene,protein_coding',
            ])

            run.main([
                os.path.join(self.test_data_directory, "input.vcf"),
                'sample.name',
                'DRB1*11:01',
                'NNalign',
                output_dir.name,
                '-e2', '15',
                '--top-score-metric=lowest',
                '--top-score-metric2=ic50',
                '--keep-tmp-files',
                '-d', 'full',
                '--run-reference-proteome-similarity',
                '--peptide-fasta', self.peptide_fasta,
                '--biotypes', 'IG_V_gene,protein_coding',
                '--allow-incomplete-transcripts'
            ])
            close_mock_fhs()

            #Class I output files
            for file_name in (
                'sample.name.MHC_I.all_epitopes.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.tsv',
                'sample.name.tsv_1-24',
                'sample.name.fasta',
                'sample.name.MHC_I.all_epitopes.aggregated.tsv',
                'sample.name.MHC_I.all_epitopes.aggregated.metrics.json',
                'sample.name.MHC_I.all_epitopes.aggregated.tsv.reference_matches',
                'sample.name.MHC_I.filtered.tsv',
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
                'sample.name.MHC_II.all_epitopes.aggregated.tsv',
                'sample.name.MHC_II.all_epitopes.aggregated.metrics.json',
                'sample.name.MHC_II.all_epitopes.aggregated.tsv.reference_matches',
                'sample.name.MHC_II.filtered.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', file_name)
                expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', file_name.replace('sample.name', 'Test'))
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))


            for file_name in (
                'sample.name.MHC_II.all_epitopes.tsv',
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
                    '--top-score-metric2=ic50',
                    '--keep-tmp-files',
                    '--run-reference-proteome-similarity',
                    '--peptide-fasta', self.peptide_fasta,
                    '--biotypes', 'IG_V_gene,protein_coding',
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
            '--top-score-metric', 'lowest',
            '--top-score-metric2', 'ic50',
            '--biotypes', 'IG_V_gene,protein_coding',
        ]
        run.main(params)
        output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'Test.MHC_II.filtered.tsv')
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
            '--biotypes', 'IG_V_gene,protein_coding',
            '--allow-incomplete-transcripts',
        ]
        run.main(params)

        for file_name in ['Test.8.fa.split_1-822', 'Test.8.fa.split_1-822.key']:
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
            expected_file = os.path.join(self.test_data_directory, 'phased', 'MHC_Class_I', 'tmp', file_name)
            self.assertTrue(cmp(output_file, expected_file, False), "files don't match %s - %s" %(output_file, expected_file))
        for file_name in ['Test.MHC_I.all_epitopes.tsv', 'Test.fasta']:
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
            expected_file = os.path.join(self.test_data_directory, 'phased', 'MHC_Class_I', file_name)
            self.assertTrue(compare(output_file, expected_file))
        for file_name in ['Test.proximal_variants.tsv', 'Test.MHC_I.all_epitopes.aggregated.tsv', 'Test.MHC_I.filtered.tsv']:
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
            '--top-score-metric2=ic50',
            '--keep-tmp-files',
            '-d', 'full',
            '--biotypes', 'IG_V_gene,protein_coding',
        ])

        for file_name in (
            'Test.Combined.all_epitopes.tsv',
            'Test.Combined.filtered.tsv',
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

    def test_problematic_amino_acids(self):
        output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request:
            run.main([
                os.path.join(self.test_data_directory, "input.vcf"),
                'Test',
                'HLA-G*01:09,HLA-E*01:01',
                'NetMHC',
                'PickPocket',
                output_dir.name,
                '-e1', '9,10',
                '--biotypes', 'IG_V_gene,protein_coding',
                '--problematic-amino-acids', 'C',
                '--allow-incomplete-transcripts'
            ])

        for file_name in (
            'Test.MHC_I.all_epitopes.tsv',
            'Test.MHC_I.filtered.tsv',
            'Test.MHC_I.all_epitopes.aggregated.tsv',
        ):
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
            expected_file = os.path.join(self.test_data_directory, 'problematic_amino_acids', 'MHC_Class_I', file_name)
            self.assertTrue(compare(output_file, expected_file))
        output_dir.cleanup()

    def test_pvacseq_run_with_ml_predictions(self):
        with patch('pvactools.lib.call_iedb.requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request, patch('pvactools.lib.net_chop.NetChop.post_query', unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            'net_chop.html'
        ))), patch('pvactools.lib.netmhc_stab.NetMHCStab.query_netmhcstabpan_server',  unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            'Netmhcstab.html'
        ))):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

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
                '--top-score-metric2=ic50',
                '--keep-tmp-files',
                '--net-chop-method', 'cterm',
                '--netmhc-stab',
                '--tdna-vaf', '0.2',
                '-d', 'full',
                '--pass-only',
                '--run-reference-proteome-similarity',
                '--peptide-fasta', self.peptide_fasta,
                '--biotypes', 'IG_V_gene,protein_coding',
                '--run-ml-predictions'
            ])
            close_mock_fhs()

        # Check that ML prediction output file exists
        ml_output_file = os.path.join(output_dir.name, 'ml_predict', 'Test.MHC_I.all_epitopes.aggregated.ML_predicted.tsv')
        
        # Check if all required pvacview files exist in ml_predict directory
        pvacview_files = [
            'Test.MHC_I.all_epitopes.aggregated.ML_predicted.tsv',
            'Test.MHC_I.all_epitopes.aggregated.metrics.json',
            'Test.MHC_II.all_epitopes.aggregated.tsv',
        ]
        for file in pvacview_files:
            file_path = os.path.join(output_dir.name, 'ml_predict', file)
            self.assertTrue(os.path.exists(file_path), f"pvacview file not found: {file_path}")

        # Check that the ML prediction output file matches the expected file
        expected_file = os.path.join(self.test_data_directory, 'ml_predictor', 'Test.MHC_I.all_epitopes.aggregated.ML_predicted.tsv')
        self.assertTrue(compare(ml_output_file, expected_file))
        
        output_dir.cleanup()