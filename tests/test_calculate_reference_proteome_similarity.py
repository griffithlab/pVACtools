import unittest
import unittest.mock
import os
import tempfile
from filecmp import cmp
import sys
import py_compile
from Bio.Blast import NCBIWWW

from pvactools.lib.calculate_reference_proteome_similarity import CalculateReferenceProteomeSimilarity
from .test_utils import *

class CalculateReferenceProteomeSimilarityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "lib", "calculate_reference_proteome_similarity.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "calculate_reference_proteome_similarity")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_calculate_self_similarity(self):
        with unittest.mock.patch('Bio.Blast.NCBIWWW.qblast', side_effect=mock_ncbiwww_qblast):
            input_file = os.path.join(self.test_data_dir, 'input.tsv')
            input_fasta = os.path.join(self.test_data_dir, 'input.fasta')
            output_file = tempfile.NamedTemporaryFile()
            metric_file = "{}.reference_matches".format(output_file.name)
            self.assertFalse(CalculateReferenceProteomeSimilarity(input_file, input_fasta, output_file.name).execute())
            self.assertTrue(cmp(
                output_file.name,
                os.path.join(self.test_data_dir, "output.tsv"),
            ))
            self.assertTrue(cmp(
                metric_file,
                os.path.join(self.test_data_dir, "output.tsv.reference_matches"),
            ))
            os.remove(metric_file)
            close_mock_fhs()

    def test_blastp_db_incompatible_with_species(self):
        with self.assertRaises(Exception) as context:
            input_file = os.path.join(self.test_data_dir, 'input.tsv')
            input_fasta = os.path.join(self.test_data_dir, 'input.fasta')
            output_file = tempfile.NamedTemporaryFile()
            CalculateReferenceProteomeSimilarity(input_file, input_fasta, output_file.name, blastp_db='refseq_select_prot', species='bonobo').execute()
            self.assertTrue("refseq_select_prot blastp database is only compatible with human and mouse species." in str(context.exception))
