import unittest
import os
import tempfile
from filecmp import cmp
import py_compile

from pvactools.lib.load_gtf_data import LoadGtfData
from tests.utils import *


# python -m unittest tests/test_pvacsplice_load_gtf_data.py
class LoadGtfDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # locate the bin and test_data directories
        cls.load_gtf_data_path = os.path.join(pvactools_directory(), "pvactools", "lib", "load_gtf_data.py")
        cls.test_data_path= os.path.join(pvactools_directory(), "tests", "test_data", "pvacsplice")
        cls.gtf_file_chr1 = os.path.join(cls.test_data_path, "inputs", "Homo_sapiens.GRCh38.105_chr1.sorted.filtered.gtf")
        cls.gtf_file_gencode = os.path.join(cls.test_data_path, "inputs", "Homo_sapiens.GENCODE.GRCh37.chr6_chr9.sorted.filtered.gtf")
        cls.gtf_file_build37 = os.path.join(cls.test_data_path, 'inputs', 'Homo_sapiens.GRCh37.75_chr1.sorted.filtered.gtf')

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.load_gtf_data_path))

    def test_load_gtf_data_runs_and_produces_expected_output(self):
        output_dir = tempfile.TemporaryDirectory()
        output_file = os.path.join(output_dir.name, 'sample_gtf.tsv')
        params = {
            'gtf_file'   : self.gtf_file_chr1,
            'output_file': output_file,
            'save_gtf'   : True,
        }
        LoadGtfData(**params).execute()
        expected_file = os.path.join(self.test_data_path, 'results', 'Test.gtf.tsv')
        self.assertTrue(cmp(
                    output_file,
                    expected_file,
                    False
                ), "files don't match {} - {}".format(output_file, expected_file)
            )

        output_dir.cleanup()

    def test_load_gtf_data_build37_runs_and_produces_expected_output(self):
        output_dir = tempfile.TemporaryDirectory()
        output_file = os.path.join(output_dir.name, 'sample_gtf.build37.tsv')
        params = {
            'gtf_file'   : self.gtf_file_build37,
            'output_file': output_file,
            'save_gtf'   : True,
        }
        LoadGtfData(**params).execute()
        expected_file = os.path.join(self.test_data_path, 'results', 'Test.gtf.build37.tsv')
        self.assertTrue(cmp(
                    output_file,
                    expected_file,
                    False
                ), "files don't match {} - {}".format(output_file, expected_file)
            )

        output_dir.cleanup()

    def test_load_gtf_data_GENCODE_runs_and_produces_expected_output(self):
        output_dir = tempfile.TemporaryDirectory()
        output_file = os.path.join(output_dir.name, 'sample_gtf.GENCODE.tsv')
        params = {
            'gtf_file'   : self.gtf_file_gencode,
            'output_file': output_file,
            'save_gtf'   : True,
        }
        LoadGtfData(**params).execute()
        expected_file = os.path.join(self.test_data_path, 'results', 'Test.gtf.GENCODE.tsv')
        self.assertTrue(cmp(
                    output_file,
                    expected_file,
                    False
                ), "files don't match {} - {}".format(output_file, expected_file)
            )

        output_dir.cleanup()

    def test_load_gtf_data_allow_incomplete_transcripts(self):
        output_dir = tempfile.TemporaryDirectory()
        output_file = os.path.join(output_dir.name, 'output_with_incomplete_transcripts.tsv')
        params = {
            'gtf_file'   : self.gtf_file_gencode,
            'output_file': output_file,
            'save_gtf'   : True,
            'allow_incomplete_transcripts': True
        }
        LoadGtfData(**params).execute()
        expected_file = os.path.join(self.test_data_path, 'results', 'Test.with_incomplete_transcripts.tsv')
        self.assertTrue(cmp(
                    output_file,
                    expected_file,
                    False
                ), "files don't match {} - {}".format(output_file, expected_file)
            )

        output_dir.cleanup()
