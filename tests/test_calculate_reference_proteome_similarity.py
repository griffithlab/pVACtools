import unittest
import unittest.mock
import os
import tempfile
from filecmp import cmp
import sys
import py_compile
from Bio.Blast import NCBIWWW
from urllib.request import urlopen
from shutil import copyfileobj
from tempfile import NamedTemporaryFile

from pvactools.lib.calculate_reference_proteome_similarity import CalculateReferenceProteomeSimilarity
from tests.utils import *

class CalculateReferenceProteomeSimilarityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "lib", "calculate_reference_proteome_similarity.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "calculate_reference_proteome_similarity")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_calculate_self_similarity_with_blast(self):
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

    def test_calculate_self_similarity_with_peptide_fasta(self):
        input_file = os.path.join(self.test_data_dir, 'input.tsv')
        input_fasta = os.path.join(self.test_data_dir, 'input.fasta')
        output_file = tempfile.NamedTemporaryFile()
        metric_file = "{}.reference_matches".format(output_file.name)
        url = "http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
        with urlopen(url) as fsrc, NamedTemporaryFile(delete=False) as fdst:
            copyfileobj(fsrc, fdst)
            fdst.close()
            self.assertFalse(CalculateReferenceProteomeSimilarity(input_file, input_fasta, output_file.name, peptide_fasta=fdst.name).execute())
            os.unlink(fdst.name)
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.peptide_fasta.tsv"),
        ))
        self.assertTrue(cmp(
            metric_file,
            os.path.join(self.test_data_dir, "output.peptide_fasta.tsv.reference_matches"),
        ))
        os.remove(metric_file)

    def test_blastp_db_incompatible_with_species(self):
        with self.assertRaises(Exception) as context:
            input_file = os.path.join(self.test_data_dir, 'input.tsv')
            input_fasta = os.path.join(self.test_data_dir, 'input.fasta')
            output_file = tempfile.NamedTemporaryFile()
            CalculateReferenceProteomeSimilarity(input_file, input_fasta, output_file.name, blastp_db='refseq_select_prot', species='bonobo').execute()
            self.assertTrue("refseq_select_prot blastp database is only compatible with human and mouse species." in str(context.exception))
