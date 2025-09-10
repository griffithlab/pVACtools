import unittest
import os
import tempfile
from filecmp import cmp
import sys
import py_compile

from pvactools.lib.ml_predictor import run_ml_predictions
from tests.utils import *

class MLPredictorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        #locate the bin and test_data directories
        cls.ml_predictor_path = os.path.join(pvactools_directory(), "pvactools", "lib", "ml_predictor.py")
        cls.test_data_path = os.path.join(pvactools_directory(), "tests", "test_data", "ml_predictor")
        cls.model_artifacts_path = os.path.join(pvactools_directory(), "pvactools", "supporting_files", "ml_model_artifacts")


    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.ml_predictor_path))

    def test_ml_predictions_output(self):
        # Use a fixed output directory instead of temporary one for inspection
        output_dir = os.path.join(pvactools_directory(), "tests", "test_data", "ml_predictor", "ml_predictor_test_output")
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "Test-TumorDNA.all_epitopes.aggregated.tsv")
            mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "Test-TumorDNA.all_epitopes.tsv")
            mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "Test-TumorDNA.all_epitopes.aggregated.tsv")
            
            result = run_ml_predictions(
                mhc1_agg_file,
                mhc1_all_file, 
                mhc2_agg_file,
                self.model_artifacts_path,
                output_dir,
                'test_sample',
                0.55
            )
            
            # Check that output files were created
            self.assertTrue(os.path.exists(result[0]))
            self.assertTrue(os.path.exists(result[1]))
            
            # Check that files contain data
            import pandas as pd
            df1 = pd.read_csv(result[0], sep='\t')
            df2 = pd.read_csv(result[1], sep='\t')
            
            self.assertGreater(df1.shape[0], 0)
            self.assertGreater(df2.shape[0], 0)
            self.assertIn('Evaluation', df1.columns)
            self.assertIn('Evaluation', df2.columns)
            
            # Print file paths for easy inspection
            print(f"\nTest output files saved to:")
            print(f"File 1: {result[0]}")
            print(f"File 2: {result[1]}")
            
        except Exception as e:
            # Don't clean up on error so you can inspect what was created
            print(f"Test failed, but output directory preserved: {output_dir}")
            raise
