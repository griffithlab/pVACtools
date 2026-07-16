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

from pvactools.tools.pvacfuse import run
import pvactools.tools.pvacfuse.main as pvacfuse_main
from tests.utils import *

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
        cls.peptide_fasta = os.path.join(pvactools_directory(), "tests", "test_data", "Homo_sapiens.GRCh38.pep.short.fa.gz")
        transcript_fasta = os.path.join(cls.test_data_directory, 'Homo_sapiens.GRCh38.95.cds.all.fa.gz')
        cls.unzipped_transcript_fasta = gunzip_file(transcript_fasta)

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

    def test_pvacfuse_pipeline(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory(),
            'agfusion',
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "agfusion"),
                'sample.name',
                'HLA-A*29:02',
                'NetMHC',
                output_dir.name,
                self.unzipped_transcript_fasta,
                '-e1', '9',
                '--top-score-metric=lowest',
                '--top-score-metric2=ic50',
                '--keep-tmp-files',
                '--run-reference-proteome-similarity',
                '--peptide-fasta', self.peptide_fasta,
                '--fasta-size', '600',
            ])
            close_mock_fhs()

            for file_name in (
                'sample.name.9.fa',
                'sample.name.transcripts.fa',
                'sample.name.tsv',
            ):
                output_file   = os.path.join(output_dir.name, file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'inputs.yml',
            ):
                output_file   = os.path.join(output_dir.name, 'log', file_name)
                self.assertTrue(os.path.exists(output_file))

            for file_name in (
                'sample.name.fasta',
                'sample.name.MHC_I.all_epitopes.tsv',
                'sample.name.MHC_I.filtered.tsv',
                'sample.name.MHC_I.all_epitopes.aggregated.tsv',
                'sample.name.MHC_I.all_epitopes.aggregated.metrics.json',
                'sample.name.MHC_I.all_epitopes.aggregated.tsv.reference_matches',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file),  "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.HLA-A*29:02.9.parsed.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', '9', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', '9', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            output_dir.cleanup()

    def test_pvacfuse_pipeline_agfusion_starfusion(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory(),
            'agfusion_starfusion',
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "agfusion_HCC1395"),
                'sample.name',
                'HLA-A*29:02',
                'NetMHC',
                output_dir.name,
                self.unzipped_transcript_fasta,
                '-e1', '9',
                '--keep-tmp-files',
                '--starfusion-file', os.path.join(self.test_data_directory, 'star-fusion.fusion_predictions.abridged.tsv'),
                '--fasta-size', '600',
            ])
            close_mock_fhs()

            for file_name in (
                'sample.name.9.fa',
                'sample.name.transcripts.fa',
                'sample.name.tsv',
            ):
                output_file   = os.path.join(output_dir.name, file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions_agfusion_starfusion', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'inputs.yml',
            ):
                output_file   = os.path.join(output_dir.name, 'log', file_name)
                self.assertTrue(os.path.exists(output_file))

            for file_name in (
                'sample.name.fasta',
                'sample.name.MHC_I.all_epitopes.tsv',
                'sample.name.MHC_I.filtered.tsv',
                'sample.name.MHC_I.all_epitopes.aggregated.tsv',
                'sample.name.MHC_I.all_epitopes.aggregated.metrics.json',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions_agfusion_starfusion', 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file),  "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.HLA-A*29:02.9.parsed.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', '9', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions_agfusion_starfusion', 'MHC_Class_I', '9', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            output_dir.cleanup()

    def test_pvacfuse_pipeline_arriba(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory(),
            'arriba',
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "arriba_fusions.tsv"),
                'sample.name',
                'HLA-A*29:02',
                'NetMHC',
                output_dir.name,
                self.unzipped_transcript_fasta,
                '-e1', '9',
                '--fasta-size', '2000',
            ])
            close_mock_fhs()

            for file_name in (
                'sample.name.9.fa',
                'sample.name.transcripts.fa',
                'sample.name.tsv',
            ):
                output_file   = os.path.join(output_dir.name, file_name)
                expected_file = os.path.join(self.test_data_directory, 'arriba_fusions', file_name)
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'inputs.yml',
            ):
                output_file   = os.path.join(output_dir.name, 'log', file_name)
                self.assertTrue(os.path.exists(output_file))

            for file_name in (
                'sample.name.fasta',
                'sample.name.MHC_I.all_epitopes.tsv',
                'sample.name.MHC_I.filtered.tsv',
                'sample.name.MHC_I.all_epitopes.aggregated.tsv',
                'sample.name.MHC_I.all_epitopes.aggregated.metrics.json',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'arriba_fusions', 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file),  "files don't match %s - %s" %(output_file, expected_file))

            output_dir.cleanup()
