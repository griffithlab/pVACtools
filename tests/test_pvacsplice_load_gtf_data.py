import unittest
import os
import tempfile
from filecmp import cmp
import py_compile

from pvactools.lib.load_gtf_data import LoadGtfData
from tests.utils import *

#python -m unittest tests/test_pvacsplice_load_gtf_data.py
class LoadGtfDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.load_gtf_data_path = os.path.join(pvactools_directory(), "pvactools", "lib", "load_gtf_data.py")
        cls.test_data_path= os.path.join(pvactools_directory(), "tests", "test_data", "pvacsplice")
        cls.gtf_file = os.path.join(cls.test_data_path, 'inputs', 'Homo_sapiens.GRCh38.105_chr1.sorted.gtf.gz')

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.load_gtf_data_path))

    def test_load_gtf_data_runs_and_produces_expected_output(self):
        output_dir = tempfile.TemporaryDirectory()
        for tsl in [1, 3, 5]:
            output_file = os.path.join(output_dir.name, f'sample_gtf_{tsl}.tsv')
            params = {
                'gtf_file'   : self.gtf_file,
                'output_file': output_file,
                'save_gtf'   : True,
                'tsl'        : tsl
            }
            LoadGtfData(**params).execute
            expected_file = os.path.join(self.test_data_path, 'results', f'Test.gtf_{tsl}.tsv')
            
            self.assertTrue(cmp(
                        output_file,
                        expected_file,
                        False
                    ), "files don't match {} - {}".format(output_file, expected_file)
                )

        output_dir.cleanup()