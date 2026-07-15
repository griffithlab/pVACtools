import unittest
import os
import py_compile

from pvactools.lib.run_argument_utils import *
from tests.utils import *

#python -m unittest tests/test_run_utils.py
class RunUtilsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin
        cls.utils_path = os.path.join(pvactools_directory(), "pvactools", "lib", "run_argument_utils.py")

    def test_module_compiles(self):
        self.assertTrue(py_compile.compile(self.utils_path))

    def test_pvacsplice_anchors_checker(self):
        checker = pvacsplice_anchors()
        self.assertEqual(
            checker("A,D,NDA"),
            ["A", "D", "NDA"]
        )

        with self.assertRaises(Exception) as context:
            checker("Test,A")

        self.assertEqual("List element must be one of 'A', 'D', 'NDA', not Test", str(context.exception))

    def test_float_range_checker(self):
        checker = float_range(0.0, 100.0)
        self.assertEqual(
            checker("0.5"),
            0.5
        )

        with self.assertRaises(Exception) as context:
            checker("Test")

        self.assertEqual("must be a floating point number", str(context.exception))

        with self.assertRaises(Exception) as context:
            checker("102.0")

        self.assertEqual("must be in range [0.0 .. 100.0]", str(context.exception))
