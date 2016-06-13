import unittest
import os
import tempfile
from filecmp import cmp
from subprocess import call
import sys
import py_compile

#python -m unittest tests/test_binding_filter.py
class BindingFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.pVac_directory = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.binding_filter_path = os.path.join(cls.pVac_directory, "pvacseq", "lib", "binding_filter.py")
        cls.test_data_path= os.path.join(cls.pVac_directory, "tests", "test_data", "binding_filter")

        #locate fof files in test data directory
        fof = os.path.join(cls.test_data_path, 'Test.fof')

        #now edit the paths contained in the fof file to point to valid files
        reader   = open(fof, mode = 'r')
        temp_fof = tempfile.NamedTemporaryFile().name
        writer   = open(temp_fof, mode = 'w')
        intake   = reader.readline().rstrip()
        while intake != "":
            #if the path is already valid, just keep it
            if os.path.isfile(intake):
                writer.write(intake + "\n")
                intake = reader.readline().rstrip()
                continue
            #if not, try changing the path to the test_data_path
            intake = os.path.join(cls.test_data_path,
                                  os.path.basename(intake))
            if os.path.exists(intake):
                writer.write(intake + "\n")
            intake = reader.readline().rstrip()
        writer.close()
        reader.close()
        cls.fof = temp_fof

    def test_binding_filter_runs_and_produces_expected_output(self):
        compiled_script_path = py_compile.compile(self.binding_filter_path)
        self.assertTrue(compiled_script_path)
        output_file = tempfile.NamedTemporaryFile().name
        binding_filter_cmd = "%s  %s  %s  %s  %s" % (
            sys.executable,
            compiled_script_path,
            os.path.join(self.test_data_path, "annotated_variants.tsv"),
            self.fof,
            output_file,
            )
        self.assertFalse(call([binding_filter_cmd], shell=True))
        self.assertTrue(cmp(
            output_file,
            os.path.join(self.test_data_path, "Test_filtered.xls"),
            False
        ))
