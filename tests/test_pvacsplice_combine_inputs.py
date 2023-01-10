import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

from pvactools.lib.combine_inputs import CombineInputs
from tests.utils import *

class CombineInputsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python        = sys.executable
        cls.executable    = os.path.join(pvactools_directory(), "pvactools", "tools", "pvacsplice", "combine_inputs.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "pvacsplice", "combine_inputs")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_combine_inputs_runs_and_produces_expected_output(self):
        junctions_file = os.path.join(self.test_data_dir, 'Test_filtered.tsv')
        variant_file = os.path.join(self.test_data_dir, 'Test_annotated.tsv')
        output_dir = tempfile.TemporaryDirectory()

        # iterate free params
        for tsl in [1, 3, 5]:
            output_file   = os.path.join(output_dir.name, 'sample_{}_combined.tsv'.format(tsl))
            #output_file   = os.path.join(self.test_data_dir, 'Test_{}_combined.tsv'.format(tsl))
            params = {
                'junctions_file'                    : junctions_file,
                'variant_file'                      : variant_file,
                #'sample_name'                       : 'sample_{}'.format(tsl),
                'output_file'                       : output_file,
                'maximum_transcript_support_level'  : tsl,
            }
            combined = CombineInputs(**params)
            combined.execute()

            #direct file comparison
            for file_name in (
                'sample_{}_combined.tsv'.format(tsl),
                #'sample_{}_annotated_filtered.tsv'.format(tsl),
            ):
                output_file   = os.path.join(output_dir.name, file_name)
                expected_file = os.path.join(self.test_data_dir, file_name.replace('sample', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match {} - {}".format(output_file, expected_file))

        output_dir.cleanup()
