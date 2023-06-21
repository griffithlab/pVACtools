import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile

from pvactools.lib.identify_problematic_amino_acids import IdentifyProblematicAminoAcids
from tests.utils import *

#python -m unittest tests/test_identify_problematic_amino_acids.py
class IdentifyProblematicAminoAcidsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.identify_problematic_amino_acids_path = os.path.join(pvactools_directory(), "pvactools", "lib", "identify_problematic_amino_acids.py")
        cls.test_data_path= os.path.join(pvactools_directory(), "tests", "test_data", "identify_problematic_amino_acids")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.identify_problematic_amino_acids_path))

    def test_invalid_problematic_amino_acid_entry(self):
        with self.assertRaises(Exception) as context:
            IdentifyProblematicAminoAcids.process_problematic_aa_entry("AX")
        self.assertEqual(
            str(context.exception),
            'Error with problematic amino acid input entry "AX". Amino acid X not supported.'
        )

        with self.assertRaises(Exception) as context:
            IdentifyProblematicAminoAcids.process_problematic_aa_entry("AX:1")
        self.assertEqual(
            str(context.exception),
            'Error with problematic amino acid input entry "AX:1". Amino acid X not supported.'
        )

        with self.assertRaises(Exception) as context:
            IdentifyProblematicAminoAcids.process_problematic_aa_entry("A:no_number")
        self.assertEqual(
            str(context.exception),
            'Error with problematic amino acid input entry "A:no_number". Position no_number isn\'t an integer.',
        )

        with self.assertRaises(Exception) as context:
            IdentifyProblematicAminoAcids.process_problematic_aa_entry("A:0")
        self.assertEqual(
            str(context.exception),
            'Error with problematic amino acid input entry "A:0". Position can\'t be 0.',
        )

    def test_identify_problematic_amino_acids_runs_and_produces_expected_output(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(IdentifyProblematicAminoAcids(
            os.path.join(
                self.test_data_path,
                'Test.all_epitopes.tsv'
            ),
            output_file.name,
            ["C", "P:1", "G:-3"],
            file_type = 'pVACseq',
            filter_type = 'soft',
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.all_epitopes.problematic.tsv"),
            False
        ))

    def test_identify_problematic_amino_acids_pvacbind_runs_and_produces_expected_output(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(IdentifyProblematicAminoAcids(
            os.path.join(
                self.test_data_path,
                'Test.all_epitopes.pvacbind.tsv'
            ),
            output_file.name,
            ["C", "P:1", "G:-3"],
            file_type = 'pVACbind',
            filter_type = 'soft',
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.all_epitopes.problematic.pvacbind.tsv"),
            False
        ))

    def test_identify_problematic_amino_acids_hard_filter_runs_and_produces_expected_output(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(IdentifyProblematicAminoAcids(
            os.path.join(
                self.test_data_path,
                'Test.all_epitopes.tsv'
            ),
            output_file.name,
            ["C", "P:1", "G:-3"],
            file_type = 'pVACseq',
            filter_type = 'hard',
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.all_epitopes.problematic.hard.tsv"),
            False
        ))
