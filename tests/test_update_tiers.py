import unittest
import unittest.mock
import os
import tempfile
from filecmp import cmp
import sys
import shutil
from tempfile import NamedTemporaryFile

from pvactools.lib.update_tiers import PvacseqUpdateTiers, PvacspliceUpdateTiers, PvacbindUpdateTiers, PvacfuseUpdateTiers
from tests.utils import *

class UpdateTiersTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "lib", "update_tiers.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "update_tiers")

    def module_compiles(self):
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
        self.assertFalse(PvacspliceUpdateTiers(
            tmp_input_file.name,
            0.5
        ).execute())
        self.assertTrue(cmp(
            tmp_input_file.name,
            os.path.join(self.test_data_dir, "HCC1395.pvacsplice.all_epitopes.aggregated.out.tsv"),
        ))
        tmp_input_file.close()

    def test_update_tiers_pvacbind(self):
        input_file = os.path.join(self.test_data_dir, 'pvacbind.aggregated.tsv')
        tmp_input_file = tempfile.NamedTemporaryFile()
        shutil.copy(input_file, tmp_input_file.name)
        self.assertFalse(PvacbindUpdateTiers(
            tmp_input_file.name,
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
        self.assertFalse(PvacfuseUpdateTiers(
            tmp_input_file.name,
        ).execute())
        self.assertTrue(cmp(
            tmp_input_file.name,
            os.path.join(self.test_data_dir, "pvacfuse.aggregated.out.tsv"),
        ))
        tmp_input_file.close()
