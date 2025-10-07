import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile

from pvactools.lib.binding_filter import BindingFilter
from tests.utils import *

#python -m unittest tests/test_binding_filter.py
class BindingFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.binding_filter_path = os.path.join(pvactools_directory(), "pvactools", "lib", "binding_filter.py")
        cls.test_data_path= os.path.join(pvactools_directory(), "tests", "test_data", "binding_filter")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.binding_filter_path))

    def test_binding_filter_runs_and_produces_expected_output_top_score_metric_median_percentile_strategy_conservative(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(BindingFilter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.binding.tsv"),
            False
        ))

    def test_binding_filter_runs_and_produces_expected_output_top_score_metric_lowest_percentile_strategy_conservative(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(BindingFilter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            top_score_metric='lowest',
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.binding.tsv"),
            False
        ))

    def test_binding_filter_runs_and_produces_expected_output_top_score_metric_median_percentile_strategy_exploratory(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(BindingFilter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            percentile_threshold_strategy='exploratory',
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.binding.exploratory.tsv"),
            False
        ))

    def test_binding_filter_runs_and_produces_expected_output_top_score_metric_lowest_percentile_strategy_explorative(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(BindingFilter(
            os.path.join(
                self.test_data_path,
                'Test.combined.parsed.tsv'
            ),
            output_file.name,
            top_score_metric='lowest',
            percentile_threshold_strategy='exploratory',
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_path, "Test.filtered.binding.exploratory.tsv"),
            False
        ))
