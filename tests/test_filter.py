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

    def module_compiles(self):
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
                exclude_nas=False
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
                exclude_nas=False
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
                exclude_nas=False
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.eq.tsv"),
            False
        ))

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
                exclude_nas=False
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
                exclude_nas=False
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
                exclude_nas=False
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.NA.tsv"),
            False
        ))

    def test_exclude_NA(self):
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
                exclude_nas=True
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.exclude_NA.tsv"),
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
                exclude_nas=True
            )],
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "output.inf.tsv"),
            False
        ))
    
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
                exclude_nas=False
            ), FilterCriterion(
                "Corresponding Fold Change",
                "<",
                "16000",
                exclude_nas=False
            )],
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
                exclude_nas=False
            ), FilterCriterion(
                "Corresponding Fold Change",
                ">",
                "16000",
                exclude_nas=False
            )],
            [],
            1
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.lt.tsv"),
            False
        ))
