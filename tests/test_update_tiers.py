import unittest
import unittest.mock
import os
import csv
import tempfile
from filecmp import cmp
import sys
import shutil
from tempfile import NamedTemporaryFile
import py_compile

from pvactools.lib.update_tiers import PvacseqUpdateTiers, PvacspliceUpdateTiers, PvacbindUpdateTiers, PvacfuseUpdateTiers
from tests.utils import *

class UpdateTiersTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "lib", "update_tiers.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "update_tiers")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_update_tiers_pvacseq(self):
        input_file = os.path.join(self.test_data_dir, 'HCC1395.all_epitopes.aggregated.tsv')
        tmp_input_file = tempfile.NamedTemporaryFile()
        shutil.copy(input_file, tmp_input_file.name)
        input_metrics_file = os.path.join(self.test_data_dir, 'HCC1395.all_epitopes.aggregated.metrics.json')
        tmp_input_metrics_file = tempfile.NamedTemporaryFile()
        shutil.copy(input_metrics_file, tmp_input_metrics_file.name)
        self.assertFalse(PvacseqUpdateTiers(
            tmp_input_file.name,
            0.5,
            percentile_threshold_strategy='exploratory',
            metrics_file=tmp_input_metrics_file.name
        ).execute())
        self.assertTrue(cmp(
            tmp_input_file.name,
            os.path.join(self.test_data_dir, "HCC1395.all_epitopes.aggregated.out.tsv"),
        ))
        self.assertTrue(cmp(
            tmp_input_metrics_file.name,
            os.path.join(self.test_data_dir, "HCC1395.all_epitopes.aggregated.metrics.out.json"),
        ))
        tmp_input_file.close()
        tmp_input_metrics_file.close()

    def test_update_tiers_pvacsplice(self):
        input_file = os.path.join(self.test_data_dir, 'HCC1395.pvacsplice.all_epitopes.aggregated.tsv')
        tmp_input_file = tempfile.NamedTemporaryFile()
        shutil.copy(input_file, tmp_input_file.name)
        input_metrics_file = os.path.join(self.test_data_dir, 'HCC1395.pvacsplice.all_epitopes.aggregated.metrics.json')
        tmp_input_metrics_file = tempfile.NamedTemporaryFile()
        shutil.copy(input_metrics_file, tmp_input_metrics_file.name)
        self.assertFalse(PvacspliceUpdateTiers(
            tmp_input_file.name,
            0.5,
            metrics_file=tmp_input_metrics_file.name
        ).execute())
        self.assertTrue(cmp(
            tmp_input_file.name,
            os.path.join(self.test_data_dir, "HCC1395.pvacsplice.all_epitopes.aggregated.out.tsv"),
        ))
        self.assertTrue(cmp(
            tmp_input_metrics_file.name,
            os.path.join(self.test_data_dir, "HCC1395.pvacsplice.all_epitopes.aggregated.out.metrics.json"),
        ))
        tmp_input_file.close()
        tmp_input_metrics_file.close()

    def test_pvacseq_missing_rna_depth_is_not_low_expression(self):
        input_file = os.path.join(self.test_data_dir, 'HCC1395.all_epitopes.aggregated.tsv')
        with open(input_file) as input_fh:
            mutation = next(csv.DictReader(input_fh, delimiter='\t'))
        mutation.update({
            'IC50 MT': '1',
            'IC50 %ile MT': '1',
            'IM %ile MT': '1',
            'Pres %ile MT': '1',
            'RNA Expr': '0',
            'RNA VAF': '0.5',
            'Allele Expr': '0',
            'RNA Depth': 'NA',
            'DNA VAF': '1',
        })

        updater = PvacseqUpdateTiers(input_file, 0.5)
        updater.anchor_calculator.is_anchor_residue_pass = unittest.mock.Mock(return_value=True)
        self.addCleanup(updater.output_file.close)

        self.assertEqual(updater.get_tier(mutation), 'NoExpr')

    def test_pvacsplice_missing_rna_depth_is_not_low_expression(self):
        input_file = os.path.join(self.test_data_dir, 'HCC1395.pvacsplice.all_epitopes.aggregated.tsv')
        with open(input_file) as input_fh:
            mutation = next(csv.DictReader(input_fh, delimiter='\t'))
        mutation.update({
            'IC50 MT': '1',
            'IC50 %ile MT': '1',
            'IM %ile MT': '1',
            'Pres %ile MT': '1',
            'RNA Expr': '0',
            'RNA VAF': '0.5',
            'Allele Expr': '0',
            'RNA Depth': 'NA',
            'DNA VAF': '1',
        })

        updater = PvacspliceUpdateTiers(input_file, 0.5)
        self.addCleanup(updater.output_file.close)

        self.assertEqual(updater.get_tier(mutation), 'NoExpr')

    def test_update_tiers_pvacbind(self):
        input_file = os.path.join(self.test_data_dir, 'pvacbind.aggregated.tsv')
        tmp_input_file = tempfile.NamedTemporaryFile()
        shutil.copy(input_file, tmp_input_file.name)
        self.assertFalse(PvacbindUpdateTiers(
            tmp_input_file.name,
            percentile_threshold_strategy='exploratory'
        ).execute())
        self.assertTrue(cmp(
            tmp_input_file.name,
            os.path.join(self.test_data_dir, "pvacbind.aggregated.out.tsv"),
        ))
        tmp_input_file.close()

    def test_update_tiers_pvacfuse(self):
        input_file = os.path.join(self.test_data_dir, 'pvacfuse.aggregated.tsv')
        tmp_input_file = tempfile.NamedTemporaryFile()
        shutil.copy(input_file, tmp_input_file.name)
        input_metrics_file = os.path.join(self.test_data_dir, 'pvacfuse.aggregated.metrics.json')
        tmp_input_metrics_file = tempfile.NamedTemporaryFile()
        shutil.copy(input_metrics_file, tmp_input_metrics_file.name)
        self.assertFalse(PvacfuseUpdateTiers(
            tmp_input_file.name,
            binding_threshold=100,
            metrics_file=tmp_input_metrics_file.name
        ).execute())
        self.assertTrue(cmp(
            tmp_input_file.name,
            os.path.join(self.test_data_dir, "pvacfuse.aggregated.out.tsv"),
        ))
        self.assertTrue(cmp(
            tmp_input_metrics_file.name,
            os.path.join(self.test_data_dir, "pvacfuse.aggregated.out.metrics.json"),
        ))
        tmp_input_file.close()
        tmp_input_metrics_file.close()
