import unittest
#from unittest.mock import patch
import os
import sys
import tempfile
from filecmp import cmp
import py_compile

from pvactools.lib.create_personal_fasta import create_personal_fasta
from tests.utils import *

#python -m unittest tests/test_pvacsplice_create_personal_fasta.py
class CreatePersonalFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python = sys.executable
        cls.executable = os.path.join(pvactools_directory(), "pvactools", "lib", "create_personal_fasta.py")
        cls.test_data_dir= os.path.join(pvactools_directory(), "tests", "test_data", "pvacsplice")
        # inputs
        cls.fasta_path = os.path.join(cls.test_data_dir, 'inputs', 'all_sequences_chr1.fa.gz')
        cls.annotated_vcf = os.path.join(cls.test_data_dir, 'inputs', 'annotated.expression_chr1.vcf.gz')

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_create_personal_fasta_runs_and_produces_expected_output(self):
        output_dir = tempfile.TemporaryDirectory()
        output_file = os.path.join(output_dir.name, 'sample_alt.fa.gz')
        expected_file = os.path.join(self.test_data_dir, 'results', 'Test.all_sequences_chr1_alt.fa.gz')
        
        # test building from scratch
        create_personal_fasta(self.fasta_path, output_file, self.annotated_vcf, 'HCC1395_TUMOR_DNA')
        self.assertTrue(cmp(
                    output_file,
                    expected_file,
                    False
                ), "files don't match {} - {}".format(output_file, expected_file)
            )

        # test when file is complete and tempdir still exists
        # test by mock_print statements in test code?
        create_personal_fasta(self.fasta_path, output_file, self.annotated_vcf, 'HCC1395_TUMOR_DNA')
        self.assertTrue(os.path.exists(output_file))

        output_dir.cleanup()