import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile

from pvactools.lib.mark_genes_of_interest import MarkGenesOfInterest
from tests.utils import *

#python -m unittest tests/test_identify_problematic_amino_acids.py
class IdentifyProblematicAminoAcidsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.mark_genes_of_interest_path = os.path.join(pvactools_directory(), "pvactools", "lib", "mark_genes_of_interest.py")
        cls.test_data_path= os.path.join(pvactools_directory(), "tests", "test_data", "mark_genes_of_interest")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.mark_genes_of_interest))

    def test_mark_genes_of_interest_output(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(MarkGenesOfInterest(
            os.path.join(
                self.test_data_path,
                'Test.all_epitopes.tsv'
            ),
            output_file.name,
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.all_epitopes.genes_of_interest.tsv"),
            False
        ))

    def test_mark_genes_of_interest_custom_fileoutput(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(MarkGenesOfInterest(
            os.path.join(
                self.test_data_path,
                'Test.all_epitopes.tsv'
            ),
            output_file.name,
            os.path.join(
                self.test_data_path,
                'cancer_census_hotspot_gene_list.tsv'
            ),
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.all_epitopes.genes_of_interest.tsv"),
            False
        ))
