import sys
import os
import unittest
import tempfile
from filecmp import cmp
import py_compile
import itertools

from pvactools.lib.fusion_to_fasta import FusionToFasta
from tests.utils import *

class FastaGeneratorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'fusion_to_fasta.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'fusion_to_fasta')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_arriba_generates_expected_file(self):
        input_file = os.path.join(self.test_data_dir, 'agfusion.tsv')
        transcript_fasta = os.path.join(self.test_data_dir, 'Homo_sapiens.GRCh38.95.cds.all.fa.gz')
        unzipped_transcript_fasta = gunzip_file(transcript_fasta)
        output_file = tempfile.NamedTemporaryFile()

        generate_fasta_params = {
            'input_file'       : input_file,
            'transcript_fasta' : unzipped_transcript_fasta,
            'output_file'      : output_file.name,
        }
        generator = FusionToFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_agfusion.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))
