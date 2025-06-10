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

        # valid human allele and epitope length
        self.assertEqual(
            get_anchor_positions("HLA-A*01:01", 8, agg_obj.allele_specific_anchors, agg_obj.anchor_probabilities, agg_obj.anchor_contribution_threshold, agg_obj.anchor_probabilities),
            [8, 2, 3]
        )

    def test_pvacsplice_anchors_checker(self):
        checker = pvacsplice_anchors()
        self.assertEqual(
            checker("A,D,NDA"),
            ["A", "D", "NDA"]
        )

        with self.assertRaises(Exception) as context:
            checker("Test,A")

        self.assertEqual("List element must be one of 'A', 'D', 'NDA', 'DA', 'N', not Test", str(context.exception))

    def test_float_range_checker(self):
        checker = float_range(0.0, 100.0)
        self.assertEqual(
            checker("0.5"),
            0.5
        )

        with self.assertRaises(Exception) as context:
            checker("Test")

        self.assertEqual("must be a floating point number", str(context.exception))

        with self.assertRaises(Exception) as context:
            checker("102.0")

        self.assertEqual("must be in range [0.0 .. 100.0]", str(context.exception))
