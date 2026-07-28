import unittest
import os
import tempfile
from filecmp import cmp
import py_compile

from pvactools.lib.filter import Filter, FilterCriterion
from tests.utils import *

class FilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.filter_path = os.path.join(pvactools_directory(), "pvactools", "lib", "filter.py")
        cls.test_data_path= os.path.join(pvactools_directory(), "tests", "test_data", "filter")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.filter_path))

    def test_less_than(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(Filter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            [FilterCriterion(
                "Median MT IC50 Score",
                "<",
                "500",
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.lt.tsv"),
            False
        ))

    def test_less_or_equal(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(Filter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            [FilterCriterion(
                "Median MT IC50 Score",
                "<=",
                "500",
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.le.tsv"),
            False
        ))

    def test_equal(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(Filter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            [FilterCriterion(
                "Median MT IC50 Score",
                "==",
                "500",
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.eq.tsv"),
            False
        ))

    def test_not_equal(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_file = os.path.join(tmp_dir, "input.tsv")
            output_file = os.path.join(tmp_dir, "output.tsv")
            with open(input_file, "w") as input_fh:
                input_fh.write("Score\n499\n500\n")

            self.assertFalse(Filter(
                input_file,
                output_file,
                [FilterCriterion("Score", "!=", "500")],
            ).execute())

            with open(output_file) as output_fh:
                self.assertEqual("Score\n499\n", output_fh.read())

    def test_greater_or_equal(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(Filter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            [FilterCriterion(
                "Median MT IC50 Score",
                ">=",
                "500",
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.ge.tsv"),
            False
        ))

    def test_greater_than(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(Filter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            [FilterCriterion(
                "Median MT IC50 Score",
                ">",
                "500",
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.gt.tsv"),
            False
        ))

    def test_NA(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(Filter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            [FilterCriterion(
                "Tumor RNA Depth",
                ">",
                "100",
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.NA.tsv"),
            False
        ))

    def test_inf(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(Filter(
            os.path.join(
                self.test_data_path,
                'input.inf.tsv'
            ),
            output_file.name,
            [FilterCriterion(
                "Corresponding Fold Change",
                ">",
                "100",
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "output.inf.tsv"),
            False
        ))

    def test_invalid_operator(self):
        with self.assertRaisesRegex(ValueError, "Unsupported filter operator"):
            FilterCriterion(
                "Median MT IC50 Score",
                "__import__('os').system('true')",
                "500",
            )

    def test_malicious_value_is_not_executed(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            marker_file = os.path.join(tmp_dir, "executed")
            input_file = os.path.join(tmp_dir, "input.tsv")
            output_file = os.path.join(tmp_dir, "output.tsv")
            malicious_value = "__import__('pathlib').Path({!r}).touch() or 1".format(marker_file)
            with open(input_file, "w") as input_fh:
                input_fh.write("Score\n{}\n".format(malicious_value))

            with self.assertRaisesRegex(ValueError, "Invalid numeric value"):
                Filter(
                    input_file,
                    output_file,
                    [FilterCriterion("Score", "<", "500")],
                ).execute()

            self.assertFalse(os.path.exists(marker_file))

    def test_conservative(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(Filter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            [FilterCriterion(
                "Median MT IC50 Score",
                "<",
                "500",
            ), FilterCriterion(
                "Corresponding Fold Change",
                "<",
                "16000",
            )],
            [],
            "AND"
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.lt.tsv"),
            False
        ))

    def test_exploratory(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(Filter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            [FilterCriterion(
                "Median MT IC50 Score",
                "<",
                "500",
            ), FilterCriterion(
                "Corresponding Fold Change",
                ">",
                "16000",
            )],
            [],
            "OR"
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.lt.tsv"),
            False
        ))
