import unittest
import os
from filecmp import cmp
from subprocess import call
import sys
import py_compile

#python -m unittest tests/test_binding_filter.py
class DefaultParametersTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.pVac_directory = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.binding_filter_path = os.path.join(cls.pVac_directory, "pvac_seq", "binding_filter.py")
        cls.test_data_path= os.path.join(cls.pVac_directory, "test_data")

        #locate fof files in test data directory
        cls.fof_files = [os.path.join(cls.test_data_path, f) for f in
                          os.listdir(cls.test_data_path) if(f.endswith(".fof") and
                                os.path.isfile(os.path.join(cls.test_data_path, f)))
                          ]

        #now edit the paths contained in the fof files to point to valid files
        for fof in cls.fof_files:
            reader = open(fof, mode = 'r')
            writer = open(fof + ".temp", mode = 'w')
            intake = reader.readline().rstrip()
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

    @classmethod
    def tearDownClass(cls):
        for fof in cls.fof_files:
            os.remove(fof + ".temp")
        os.remove(os.path.join(cls.test_data_path,
                               "test_binding_filter_py_default.xls"))

    def test_binding_filter_runs_and_produces_expected_output(self):
        compiled_script_path = py_compile.compile(self.binding_filter_path)
        self.assertTrue(compiled_script_path)
        binding_filter_cmd = "%s  %s  %s  %s  %s" % (
            sys.executable,
            compiled_script_path,
            os.path.join(self.test_data_path, "annotated_variants.tsv"),
            os.path.join(self.test_data_path, "Test.fof.temp"),
            os.path.join(self.test_data_path, "test_binding_filter_py_default.xls"),
            )
        self.assertFalse(call([binding_filter_cmd], shell=True))
        self.assertTrue(cmp(
            os.path.join(self.test_data_path, "test_binding_filter_py_default.xls"),
            os.path.join(self.test_data_path, "Test_filtered.xls"),
            False
        ))
