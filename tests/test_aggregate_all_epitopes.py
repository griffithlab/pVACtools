import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile

from pvactools.lib.aggregate_all_epitopes import PvacseqAggregateAllEpitopes, UnmatchedSequenceAggregateAllEpitopes
from .test_utils import *

class AggregateAllEptiopesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "lib", "aggregate_all_epitopes.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "aggregate_all_epitopes")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_aggregate_all_epitopes_pvacseq_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(PvacseqAggregateAllEpitopes(os.path.join(self.test_data_dir, 'Test.all_epitopes.tsv'), output_file.name).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertTrue(cmp(
            metrics_file,
            os.path.join(self.test_data_dir, "output.metrics.json"),
        ))
        os.remove(metrics_file)

        for i in ["ui.R", "app.R", "server.R", "styling.R", "anchor_and_helper_functions.R"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            self.assertTrue(os.path.isfile(pvacview_file))
            os.remove(pvacview_file)

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            self.assertTrue(os.path.isfile(pvacview_file))
            os.remove(pvacview_file)

    def test_aggregate_all_epitopes_pvacfuse_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(UnmatchedSequenceAggregateAllEpitopes(os.path.join(self.test_data_dir, 'Test.all_epitopes.pvacfuse.tsv'), output_file.name).execute())
        self.assertTrue(compare(
            output_file.name,
            os.path.join(self.test_data_dir, "output.pvacfuse.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertFalse(os.path.isfile(metrics_file))

        for i in ["ui.R", "app.R", "server.R", "styling.R", "anchor_and_helper_functions.R"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            self.assertFalse(os.path.isfile(pvacview_file))

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            self.assertFalse(os.path.isfile(pvacview_file))

    def test_aggregate_all_epitopes_pvacbind_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(UnmatchedSequenceAggregateAllEpitopes(os.path.join(self.test_data_dir, 'Test.all_epitopes.pvacbind.tsv'), output_file.name).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.pvacbind.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertFalse(os.path.isfile(metrics_file))

        for i in ["ui.R", "app.R", "server.R", "styling.R", "anchor_and_helper_functions.R"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            self.assertFalse(os.path.isfile(pvacview_file))

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            self.assertFalse(os.path.isfile(pvacview_file))
