import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile

from pvactools.lib.aggregate_all_epitopes import PvacseqAggregateAllEpitopes, PvacfuseAggregateAllEpitopes, PvacbindAggregateAllEpitopes, PvacspliceAggregateAllEpitopes
from tests.utils import *

class AggregateAllEpitopesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "lib", "aggregate_all_epitopes.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "aggregate_all_epitopes")
        cls.pvacview_r_files = ["ui.R", "app.R", "server.R", "styling.R", "anchor_and_helper_functions.R", "neofox_ui.R", "custom_ui.R"]

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

        for i in self.pvacview_r_files:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            self.assertTrue(os.path.isfile(pvacview_file))
            os.remove(pvacview_file)

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            self.assertTrue(os.path.isfile(pvacview_file))
            os.remove(pvacview_file)

    def test_aggregate_all_epitopes_HCC1395_pvacseq_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(PvacseqAggregateAllEpitopes(os.path.join(self.test_data_dir, 'HCC1395_TUMOR_DNA.all_epitopes.short.tsv'), output_file.name).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "HCC1395.output.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertTrue(cmp(
            metrics_file,
            os.path.join(self.test_data_dir, "HCC1395.output.metrics.json"),
        ))
        os.remove(metrics_file)

        for i in self.pvacview_r_files:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            self.assertTrue(os.path.isfile(pvacview_file))
            os.remove(pvacview_file)

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            self.assertTrue(os.path.isfile(pvacview_file))
            os.remove(pvacview_file)

    def test_aggregate_all_epitopes_pvacseq_na_mutation_position_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(PvacseqAggregateAllEpitopes(os.path.join(self.test_data_dir, 'Test.all_epitopes.na_mutation_position.tsv'), output_file.name).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.na_mutation_position.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertTrue(cmp(
            metrics_file,
            os.path.join(self.test_data_dir, "output.na_mutation_position.metrics.json"),
        ))
        os.remove(metrics_file)

        for i in self.pvacview_r_files:
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
        self.assertFalse(PvacfuseAggregateAllEpitopes(os.path.join(self.test_data_dir, 'Test.all_epitopes.pvacfuse.tsv'), output_file.name).execute())
        self.assertTrue(compare(
            output_file.name,
            os.path.join(self.test_data_dir, "output.pvacfuse.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertFalse(os.path.isfile(metrics_file))

        for i in self.pvacview_r_files:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            self.assertFalse(os.path.isfile(pvacview_file))

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            self.assertFalse(os.path.isfile(pvacview_file))

    def test_aggregate_all_epitopes_pvacbind_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(PvacbindAggregateAllEpitopes(os.path.join(self.test_data_dir, 'Test.all_epitopes.pvacbind.tsv'), output_file.name).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.pvacbind.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertFalse(os.path.isfile(metrics_file))

        for i in self.pvacview_r_files:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            self.assertFalse(os.path.isfile(pvacview_file))

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            self.assertFalse(os.path.isfile(pvacview_file))

    def test_aggregate_all_epitopes_pvacsplice_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(PvacspliceAggregateAllEpitopes(
            os.path.join(self.test_data_dir, 'Test.all_epitopes.pvacsplice.tsv'),
            output_file.name,
            aggregate_inclusion_binding_threshold=30000,
            binding_threshold=15000,
            expn_val=2.4,
        ).execute())
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.pvacsplice.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertFalse(os.path.isfile(metrics_file))

        for i in ["ui.R", "app.R", "server.R", "styling.R", "anchor_and_helper_functions.R"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            self.assertFalse(os.path.isfile(pvacview_file))

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            self.assertFalse(os.path.isfile(pvacview_file))

    def test_aggregate_all_epitopes_allele_specific_binding_thresholds_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(
            PvacseqAggregateAllEpitopes(
                os.path.join(self.test_data_dir, 'Test.all_epitopes.tsv'),
                output_file.name,
                allele_specific_binding_thresholds=True
            ).execute()
        )
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.allele_specific_binding_thresholds.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertTrue(cmp(
            metrics_file,
            os.path.join(self.test_data_dir, "output.allele_specific_binding_thresholds.metrics.json"),
        ))
        os.remove(metrics_file)

        for i in self.pvacview_r_files:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            os.remove(pvacview_file)

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            os.remove(pvacview_file)

    def test_aggregate_all_epitopes_lowest_top_score_metric_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(
            PvacseqAggregateAllEpitopes(
                os.path.join(self.test_data_dir, 'Test.all_epitopes.tsv'),
                output_file.name,
                top_score_metric="lowest",
                top_score_metric2="ic50"
            ).execute()
        )
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.lowest_top_score_metric.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertTrue(cmp(
            metrics_file,
            os.path.join(self.test_data_dir, "output.lowest_top_score_metric.metrics.json"),
        ))
        os.remove(metrics_file)

        for i in self.pvacview_r_files:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            os.remove(pvacview_file)

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            os.remove(pvacview_file)

    def test_aggregate_all_epitopes_pvacseq_top_score_metric2_percentile_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(
            PvacseqAggregateAllEpitopes(
                os.path.join(self.test_data_dir, 'Test.all_epitopes.tsv'),
                output_file.name,
                top_score_metric="lowest",
                top_score_metric2="percentile"
            ).execute()
        )
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.percentile_top_score_metric2.tsv"),
        ))

        metrics_file = output_file.name.replace('.tsv', '.metrics.json')
        self.assertTrue(cmp(
            metrics_file,
            os.path.join(self.test_data_dir, "output.percentile_top_score_metric2.metrics.json"),
        ))
        os.remove(metrics_file)

        for i in self.pvacview_r_files:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), i)
            os.remove(pvacview_file)

        for i in ["anchor.jpg", "pVACview_logo.png", "pVACview_logo_mini.png"]:
            pvacview_file = os.path.join(os.path.dirname(output_file.name), "www", i)
            os.remove(pvacview_file)

    def test_aggregate_all_epitopes_pvacbind_lowest_top_score_metric_runs_and_produces_expected_output(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        self.assertFalse(
            PvacbindAggregateAllEpitopes(
                os.path.join(self.test_data_dir, 'Test.all_epitopes.pvacbind.tsv'),
                output_file.name,
                top_score_metric="lowest",
                top_score_metric2 = "ic50", 
            ).execute()
        )
        
        self.assertTrue(cmp(
            output_file.name,
            os.path.join(self.test_data_dir, "output.pvacbind.lowest_top_score_metric.tsv"),
        ))

    def test_aggregate_all_epitopes_pvacbind_top_score_metric2_percentile(self):
        self.assertTrue(py_compile.compile(self.executable))
        output_file = tempfile.NamedTemporaryFile(suffix='.tsv')
        output_file_name = output_file.name
        self.assertFalse(
            PvacbindAggregateAllEpitopes(
                os.path.join(self.test_data_dir, "Test.all_epitopes.pvacbind.tsv"),
                output_file_name,
                top_score_metric="lowest",
                top_score_metric2="percentile"
            ).execute()
        )
        self.assertTrue(cmp(
            output_file_name,
            os.path.join(self.test_data_dir, "output.pvacbind.top_score_metric2_percentile.tsv")
        ))