import unittest
import unittest.mock
import os
import tempfile
from filecmp import cmp
import sys
import py_compile
from lib.calculate_reference_proteome_similarity import *
from Bio.Blast import NCBIWWW
from .test_utils import *

class CalculateReferenceProteomeSimilarityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.base_dir      = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.executable    = os.path.join(cls.base_dir, "lib", "calculate_reference_proteome_similarity.py")
        cls.test_data_dir = os.path.join(cls.base_dir, "tests", "test_data", "calculate_reference_proteome_similarity")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_calculate_self_similarity(self):
        with unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            input_file = os.path.join(self.test_data_dir, 'input.tsv')
            input_fasta = os.path.join(self.test_data_dir, 'input.fasta')
            output_file = tempfile.NamedTemporaryFile()
            self.assertFalse(CalculateReferenceProteomeSimilarity(input_file, input_fasta, output_file.name).execute())
            self.assertTrue(cmp(
                output_file.name,
                os.path.join(self.test_data_dir, "output.tsv"),
            ))
            self.assertTrue(cmp(
                "{}.reference_matches".format(output_file.name),
                os.path.join(self.test_data_dir, "output.tsv.reference_matches"),
            ))
            close_mock_fhs()
