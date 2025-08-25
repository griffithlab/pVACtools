import unittest
import os
import py_compile

from pvactools.lib.aggregate_all_epitopes import AnchorResiduePass
from tests.utils import *

class AnchorResiduePassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin
        cls.utils_path = os.path.join(pvactools_directory(), "pvactools", "lib", "anchor_residue_pass.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "anchor_residue_pass")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.utils_path))

    def test_get_anchor_positions(self):
        anchor_obj = AnchorResiduePass(
            binding_threshold=500,
            use_allele_specific_binding_thresholds=False,
            allele_specific_binding_thresholds=None,
            allele_specific_anchors=True,
            anchor_contribution_threshold=0.8,
        )

        # anchor positions are stored for the given valid mouse allele and epitope length
        self.assertEqual(
            anchor_obj.get_anchor_positions("H-2-Db", 11),
            [5, 11]
        )

        # valid mouse allele, but positions not stored for given epitope length
        self.assertEqual(
            anchor_obj.get_anchor_positions("H-2-Kb", 11),
            [1, 2, 10, 11]
        )

        # valid human allele and epitope length
        self.assertEqual(
            anchor_obj.get_anchor_positions("HLA-A*01:01", 8),
            [8, 2, 3]
        )
