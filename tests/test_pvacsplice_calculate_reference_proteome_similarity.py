import unittest
import os
import re
import sys
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run
from tempfile import NamedTemporaryFile

from pvactools.tools.pvacsplice import *
from tests.utils import *

def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacsplice'
    )

class PvacspliceCalculateReferenceProteomeSimilarityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pvactools_directory = pvactools_directory()
        cls.test_data_directory = test_data_directory()
        cls.peptide_fasta = os.path.join(pvactools_directory(), "tests", "test_data", "Homo_sapiens.GRCh38.pep.all.fa.gz")

    def test_command(self):
        pvac_script_path = os.path.join(
            self.pvactools_directory,
            'pvactools',
            'tools',
            'pvacsplice',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        result = subprocess_run([
            sys.executable,
            pvac_script_path,
            'calculate_reference_proteome_similarity',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacsplice calculate_reference_proteome_similarity -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            'pvactools',
            "tools",
            "pvacsplice",
            "calculate_reference_proteome_similarity.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_runs(self):
        input_file = os.path.join(self.test_data_directory, 'results', 'run', 'MHC_Class_I', 'HCC1395_TUMOR_DNA.all_epitopes.aggregated.tsv')
        input_fasta = os.path.join(self.test_data_directory, 'results', 'run', 'HCC1395_TUMOR_DNA.transcripts.fa')
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(calculate_reference_proteome_similarity.main([
            input_file,
            input_fasta,
            output_file.name,
            '--peptide-fasta', self.peptide_fasta
        ]))
