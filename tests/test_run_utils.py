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

    def test_module_compiles(self):
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

    def test_pvacsplice_anchors_checker(self):
        checker = pvacsplice_anchors()
        self.assertEqual(
            checker("A,D,NDA"),
            ["A", "D", "NDA"]
        )

        with self.assertRaises(Exception) as context:
            checker("Test,A")

        self.assertEqual("List element must be one of 'A', 'D', 'NDA', not Test", str(context.exception))

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

    def test_is_preferred_transcript_accepts_boolean_values_and_strings(self):
        for column, strategy in [
            ('Canonical', 'canonical'),
            ('MANE Select', 'mane_select'),
        ]:
            for value, expected in [
                (True, True),
                (False, False),
                ('True', True),
                ('False', False),
                ('Not Run', True),
            ]:
                with self.subTest(column=column, value=value):
                    mutation = {
                        'Canonical': False,
                        'MANE Select': False,
                    }
                    mutation[column] = value
                    self.assertEqual(
                        expected,
                        is_preferred_transcript(mutation, [strategy], 1),
                    )

    def test_is_preferred_transcript_rejects_invalid_boolean_value(self):
        for column, strategy in [
            ('Canonical', 'canonical'),
            ('MANE Select', 'mane_select'),
        ]:
            with self.subTest(column=column):
                mutation = {
                    'Canonical': False,
                    'MANE Select': False,
                }
                mutation[column] = 'yes'
                with self.assertRaises(ValueError) as context:
                    is_preferred_transcript(mutation, [strategy], 1)
                self.assertEqual(
                    "Invalid value 'yes' for {!r}. Expected True, False, or 'Not Run'.".format(column),
                    str(context.exception),
                )

    def test_is_preferred_transcript_does_not_execute_malicious_value(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            marker_file = os.path.join(tmp_dir, "executed")
            mutation = {
                'Canonical': "__import__('pathlib').Path({!r}).touch() or True".format(marker_file),
                'MANE Select': False,
            }
            with self.assertRaisesRegex(ValueError, "Invalid value"):
                is_preferred_transcript(mutation, ['canonical'], 1)
            self.assertFalse(os.path.exists(marker_file))
