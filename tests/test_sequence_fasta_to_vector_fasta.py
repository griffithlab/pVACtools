import os
import unittest
import tempfile
from filecmp import cmp
import py_compile
import itertools

from pvactools.lib.sequence_fasta_to_vector_fasta import SequenceFastaToVectorFasta
from tests.utils import *

class SequenceFastaToVectorFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'sequence_fasta_to_vector_fasta.py')
        cls.test_data_dir  = os.path.join(pvactools_directory(), 'tests', 'test_data', 'sequence_fasta_to_vector_fasta')
        cls.epitope_length = 8
        cls.flanking_sequence_length = 10

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_pvacvector_input_fasta_generates_expected_file(self):
        input_file = os.path.join(self.test_data_dir, 'pvacvector.fa')
        output_file = tempfile.NamedTemporaryFile()

        peptides = ['MT.FAT3.R4848T', 'MT.PRDM15.G654W', 'MT.DTX3L.G501R', 'MT.SUMF2.G23A', 'MT.POM121C.G3107R', 'MT.PEX1.V356I', 'MT.NRCAM.P838H', 'MT.CASP10.S654R', 'MT.ACSL3.S345N', 'MT.TP53.R157H']
        junctions_to_test = list(itertools.permutations(peptides, 2))

        generate_fasta_params = {
            'input_file'        : input_file,
            'epitope_length'    : 8,
            'output_file'       : output_file.name,
            'spacer'            : 'HH',
            'junctions_to_test' : junctions_to_test,
            'clip_length'       : 1
        }
        generator = SequenceFastaToVectorFasta(**generate_fasta_params)

        self.assertFalse(generator.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_pvacvector.fasta')
        self.assertTrue(cmp(output_file.name, expected_output_file))
