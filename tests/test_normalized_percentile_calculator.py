import unittest
import os
import py_compile
import requests

from pvactools.lib.normalized_percentile_calculator import NormalizedPercentileCalculator
from tests.utils import *

class OutputParserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        executable_dir    = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable    = os.path.join(executable_dir, 'normalized_percentile_calculator.py')
        url = f"https://raw.githubusercontent.com/griffithlab/pvactools_percentiles_data/main/hdf5/HLA-A_02_01_percentiles.h5"
        response = requests.get(url, stream=True)
        response.raise_for_status()
        cls.reference_file = "/tmp/HLA-A_02_01_percentiles.h5"
        with open(cls.reference_file, "wb") as fh:
            for chunk in response.iter_content(chunk_size=8192):
                fh.write(chunk)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.reference_file)

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_allele_normalization(self):
        calculator = NormalizedPercentileCalculator(reference_scores_path="/tmp")
        self.assertEqual(calculator.normalize_allele("HLA-A*02:01"), "HLA-A_02_01")
        self.assertEqual(calculator.normalize_allele("HLA-A*02:122"), "HLA-A_02_122")
        self.assertEqual(calculator.normalize_allele("HLA-A*02:125N"), "HLA-A_02_125N")
        self.assertEqual(calculator.normalize_allele("HLA-A*02:53N"), "HLA-A_02_53N")

    def test_normalized_percentile(self):
        calculator = NormalizedPercentileCalculator(reference_scores_path="/tmp")

        ####Per-Length Tests####
        #test specific percentile
        self.assertEqual(calculator.calculate_normalized_percentile("HLA-A_02_01", 8, 500, "NetMHC", mode="per_length"), 0.658)
        #test percentile for a different length should differ
        self.assertNotEqual(calculator.calculate_normalized_percentile("HLA-A_02_01", 8, 500, "NetMHC", mode="per_length"), calculator.calculate_normalized_percentile("HLA-A_02_01", 9, 500, "NetMHC", mode="per_length"))
        #test for unsupported length should yield "NA"
        self.assertEqual(calculator.calculate_normalized_percentile("HLA-A_02_01", 17, 500, "NetMHC", mode="per_length"), "NA")
        #test is_reversed
        self.assertEqual(calculator.calculate_normalized_percentile("HLA-A_02_01", 9, 0.9, "DeepImmuno", is_reversed=True, mode="per_length"), 2.084)

        ####Length-Agnostic Tests####
        #test specific percentile
        self.assertEqual(calculator.calculate_normalized_percentile("HLA-A_02_01", 8, 500, "NetMHC", mode="length_agnostic"), 2.348)
        #test percentile for a different length should be the same
        self.assertEqual(calculator.calculate_normalized_percentile("HLA-A_02_01", 8, 500, "NetMHC", mode="length_agnostic"), calculator.calculate_normalized_percentile("HLA-A_02_01", 9, 500, "NetMHC", mode="length_agnostic"))
        #test for unsupported length should yield a value
        self.assertEqual(calculator.calculate_normalized_percentile("HLA-A_02_01", 17, 500, "NetMHC", mode="length_agnostic"), 2.348)
        #test is_reversed
        self.assertEqual(calculator.calculate_normalized_percentile("HLA-A_02_01", 9, 0.9, "DeepImmuno", is_reversed=True, mode="length_agnostic"), 2.078)
