import unittest
import os
import py_compile

from pvactools.lib.run_utils import *
from tests.utils import *

#python -m unittest tests/test_run_utils.py
class RunUtilsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin
        cls.utils_path = os.path.join(pvactools_directory(), "pvactools", "lib", "run_utils.py")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.utils_path))
