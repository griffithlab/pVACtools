import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile

from pvactools.lib.aggregate_all_epitopes import PvacseqAggregateAllEpitopes
from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *
from tests.utils import *

#python -m unittest tests/test_run_utils.py
class RunUtilsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin
        cls.utils_path = os.path.join(pvactools_directory(), "pvactools", "lib", "run_utils.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "run_utils")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.utils_path))

    def test_combine_class_ii_alleles(self):
        self.assertEqual(
            sorted(combine_class_ii_alleles(["DQA1*06:02", "DQB1*06:43", "DQB1*06:44"])),
            sorted(["DQA1*06:02", "DQB1*06:43", "DQB1*06:44", "DQA1*06:02-DQB1*06:43", "DQA1*06:02-DQB1*06:44"])
        )

        #two different types of alleles
        self.assertEqual(
            sorted(combine_class_ii_alleles(["DQA1*06:02", "DQB1*06:43", "DPA1*01:03", "DPB1*02:01"])),
            sorted(["DQA1*06:02", "DQB1*06:43", "DPA1*01:03", "DPB1*02:01", "DQA1*06:02-DQB1*06:43", "DPA1*01:03-DPB1*02:01"])
        )

        #DRA*01:01-DRB9*01:02 is invalid
        self.assertEqual(
            sorted(combine_class_ii_alleles(["DRB9*01:02", "DRA*01:01"])),
            sorted(["DRB9*01:02", "DRA*01:01"])
        )
    
    def test_get_anchor_positions(self):
        agg_obj = PvacseqAggregateAllEpitopes(input_file=os.path.join(self.test_data_dir, 'M_GC-OxParp_A-OxParp_A_FF_DNA.all_epitopes.tsv'), 
                                                   output_file="",
                                                   allele_specific_anchors=True
                                                  )

        # anchor positions are stored for the given valid mouse allele and epitope length
        self.assertEqual(
            get_anchor_positions("H-2-Db", 11, agg_obj.allele_specific_anchors, agg_obj.anchor_probabilities, agg_obj.anchor_contribution_threshold, agg_obj.mouse_anchor_positions),
            [5, 11]
        )

        # valid mouse allele, but positions not stored for given epitope length
        self.assertEqual(
            get_anchor_positions("H-2-Kb", 11, agg_obj.allele_specific_anchors, agg_obj.anchor_probabilities, agg_obj.anchor_contribution_threshold, agg_obj.mouse_anchor_positions),
            [1, 2, 10, 11]
        )
