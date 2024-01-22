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
        cls.peptide_fasta = os.path.join(pvactools_directory(), "tests", "test_data", "Homo_sapiens.GRCh38.pep.all.fa.gz")

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
        self.assertFalse(CalculateReferenceProteomeSimilarity(input_file, input_fasta, output_file.name, peptide_fasta=self.peptide_fasta).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.peptide_fasta.tsv"),
        ))
        self.assertTrue(cmp(
            metric_file,
            os.path.join(self.test_data_dir, "output.peptide_fasta.tsv.reference_matches"),
        ))
        os.remove(metric_file)

    def test_calculate_self_similarity_with_aggregated_tsv_and_peptide_fasta(self):
        input_file = os.path.join(self.test_data_dir, 'Test.all_epitopes.aggregated.tsv')
        input_aggregated_metrics_file = os.path.join(self.test_data_dir, 'Test.all_epitopes.aggregated.tsv.metrics.json')
        tmp_aggregated_metrics_file = tempfile.NamedTemporaryFile()
        import shutil
        shutil.copy(input_aggregated_metrics_file, tmp_aggregated_metrics_file.name)
        input_fasta = os.path.join(self.test_data_dir, 'Test.fasta')
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        metric_file = "{}.reference_matches".format(output_file.name)
        output_aggregated_metrics_file = output_file.name.replace(".tsv", ".metrics.json")
        self.assertFalse(CalculateReferenceProteomeSimilarity(
            input_file,
            input_fasta,
            output_file.name,
            peptide_fasta=self.peptide_fasta,
            aggregate_metrics_file=tmp_aggregated_metrics_file.name,
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.aggregated.peptide_fasta.tsv"),
        ))
        self.assertTrue(cmp(
            metric_file,
            os.path.join(self.test_data_dir, "output.aggregated.peptide_fasta.tsv.reference_matches"),
        ))
        self.assertTrue(cmp(
            output_aggregated_metrics_file,
            os.path.join(self.test_data_dir, "output.aggregated.peptide_fasta.tsv.metrics.json"),
        ))
        os.remove(metric_file)

    def test_calculate_self_similarity_with_aggregated_tsv_and_peptide_fasta_mouse(self):
        input_file = os.path.join(self.test_data_dir, 'Test.all_epitopes.aggregated.mouse.tsv')
        input_aggregated_metrics_file = os.path.join(self.test_data_dir, 'Test.all_epitopes.aggregated.mouse.tsv.metrics.json')
        tmp_aggregated_metrics_file = tempfile.NamedTemporaryFile()
        import shutil
        shutil.copy(input_aggregated_metrics_file, tmp_aggregated_metrics_file.name)
        input_fasta = os.path.join(self.test_data_dir, 'Test.mouse.fasta')
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        metric_file = "{}.reference_matches".format(output_file.name)
        output_aggregated_metrics_file = output_file.name.replace(".tsv", ".metrics.json")
        self.assertFalse(CalculateReferenceProteomeSimilarity(
            input_file,
            input_fasta,
            output_file.name,
            peptide_fasta=self.peptide_fasta,
            aggregate_metrics_file=tmp_aggregated_metrics_file.name,
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.aggregated.peptide_fasta.mouse.tsv"),
        ))
        os.remove(metric_file)

    def test_wt_peptide_fully_in_mt_peptide(self):
        input_file = os.path.join(self.test_data_dir, 'input_wt_in_mt.tsv')
        input_fasta = os.path.join(self.test_data_dir, 'input_wt_in_mt.fasta')
        output_file = tempfile.NamedTemporaryFile()
        metric_file = "{}.reference_matches".format(output_file.name)
        self.assertFalse(CalculateReferenceProteomeSimilarity(
            input_file,
            input_fasta,
            output_file.name,
            peptide_fasta=self.peptide_fasta,
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output_wt_in_mt.tsv"),
        ))
        self.assertTrue(cmp(
            metric_file,
            os.path.join(self.test_data_dir, "output_wt_in_mt.tsv.reference_matches"),
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
