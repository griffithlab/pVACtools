import unittest
import unittest.mock
import os
import tempfile
from filecmp import cmp
import sys
import py_compile
from Bio.Blast import NCBIWWW
from urllib.request import urlopen
from shutil import copyfileobj
from tempfile import NamedTemporaryFile
import logging
from testfixtures import LogCapture, StringComparison as S

from pvactools.lib.aggregate_report_filter import AggregateReportFilter
from tests.utils import *

class AggregateReportFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "lib", "aggregate_report_filter.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "aggregate_report_filter")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_aggregate_report_filter_pvacseq(self):
        input_file = os.path.join(self.test_data_dir, 'Test.all_epitopes.aggregated.tsv')
        input_metrics_file = os.path.join(self.test_data_dir, 'Test.all_epitopes.aggregated.metrics.json')
        output_file = tempfile.NamedTemporaryFile()
        output_metrics_file = tempfile.NamedTemporaryFile()
        self.assertFalse(AggregateReportFilter(
            input_file,
            output_file.name,
            input_metrics_file,
            output_metrics_file.name,
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.aggregated.tsv"),
        ))
        self.assertTrue(cmp(
            output_metrics_file.name,
            os.path.join(self.test_data_dir, "output.aggregated.metrics.json"),
        ))

    def test_aggregate_report_filter_pvacfuse(self):
        input_file = os.path.join(self.test_data_dir, 'Test.pvacfuse.all_epitopes.aggregated.tsv')
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(AggregateReportFilter(
            input_file,
            output_file.name,
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.pvacfuse.aggregated.tsv"),
        ))
