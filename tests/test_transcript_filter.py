import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile

from pvactools.lib.transcript_filter import TranscriptFilter
from tests.utils import *

#python -m unittest tests/test_transcript_filter.py
class TranscriptFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "lib", "transcript_filter.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "transcript_filter")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_transcript_filter_runs_and_produces_expected_output(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(TranscriptFilter(
            os.path.join(
                self.test_data_dir,
                'Test.all_epitopes.tsv'
            ),
            output_file.name,
            transcript_prioritization_strategy=['tsl']
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "Test.filtered.default.tsv"),
        ))

    def test_transcript_filter_with_defaults_runs_and_produces_expected_output(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(TranscriptFilter(
            os.path.join(
                self.test_data_dir,
                'Test.all_epitopes.tsv'
            ),
            output_file.name,
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "Test.filtered.defaults.tsv"),
        ))

    def test_transcript_filter_with_max_tsl_runs_and_produces_expected_output(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(TranscriptFilter(
            os.path.join(
                self.test_data_dir,
                'Test.all_epitopes.tsv'
            ),
            output_file.name,
            maximum_transcript_support_level=3,
            transcript_prioritization_strategy=['tsl']
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "Test.filtered.max_tsl_3.tsv"),
        ))

    def test_transcript_filter_tsl_not_supported(self):
        output_file = tempfile.NamedTemporaryFile()
        self.assertFalse(TranscriptFilter(
            os.path.join(
                self.test_data_dir,
                'Test.all_epitopes.GRCh37.tsv'
            ),
            output_file.name,
            transcript_prioritization_strategy=['tsl']
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "Test.all_epitopes.GRCh37.tsv"),
        ))
