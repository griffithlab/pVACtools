import unittest
import os
import tempfile
import shutil
import pandas as pd
import py_compile
from unittest.mock import patch, MagicMock

from pvactools.lib.ml_predictor import (
    run_ml_predictions,
    merge_and_prepare_data,
    clean_and_impute_data,
    make_ml_predictions,
    create_final_output,
    _resolve_artifact_paths,
    _get_default_artifacts_dir,
    define_add_ml_predictions_parser
)
from tests.utils import *


class MLPredictorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ml_predictor_path = os.path.join(pvactools_directory(), "pvactools", "lib", "ml_predictor.py")
        cls.test_data_path = os.path.join(pvactools_directory(), "tests", "test_data", "ml_predictor")
        cls.model_artifacts_path = os.path.join(pvactools_directory(), "pvactools", "supporting_files", "ml_model_artifacts")

    def module_compiles(self):
        """Test that the module compiles without syntax errors."""
        self.assertTrue(py_compile.compile(self.ml_predictor_path))

    # ============================================================================
    # UNIT TESTS FOR INDIVIDUAL FUNCTIONS
    # ============================================================================

    def test_resolve_artifact_paths(self):
        """Test that artifact paths are resolved correctly."""
        test_dir = "/test/path"
        model_path, imputer_path, encoders_path = _resolve_artifact_paths(test_dir)
        
        self.assertEqual(model_path, "/test/path/rf_downsample_model_numpy126.pkl")
        self.assertEqual(imputer_path, "/test/path/trained_imputer_numpy126.joblib")
        self.assertEqual(encoders_path, "/test/path/label_encoders_numpy126.pkl")

    def test_merge_and_prepare_data_basic(self):
        """Test merge_and_prepare_data with valid inputs."""
        mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
        mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
        mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
        
        result = merge_and_prepare_data(mhc1_agg_file, mhc1_all_file, mhc2_agg_file)
        
        # Check result is a DataFrame
        self.assertIsInstance(result, pd.DataFrame)
        # Check required columns exist
        self.assertIn("ID", result.columns)
        self.assertIn("Evaluation", result.columns)
        self.assertIn("TSL", result.columns)
        # Check that renamed IC50/percentile columns exist (Best Peptide/Transcript are used for merging but dropped in final output)
        self.assertIn("IC50 MT class1", result.columns)
        self.assertIn("IC50 MT class2", result.columns)
        self.assertIn("%ile MT class1", result.columns)
        self.assertIn("%ile MT class2", result.columns)
        # Check data types
        self.assertTrue(pd.api.types.is_bool_dtype(result["Prob match"]))
        self.assertTrue(pd.api.types.is_bool_dtype(result["Gene of Interest"]))
        self.assertTrue(pd.api.types.is_integer_dtype(result["TSL"]))

    def test_merge_and_prepare_data_pos_transformation(self):
        """Test that Pos column transformation handles various formats."""
        mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
        mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
        mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
        
        result = merge_and_prepare_data(mhc1_agg_file, mhc1_all_file, mhc2_agg_file)
        
        # Check Pos is numeric
        self.assertTrue(pd.api.types.is_numeric_dtype(result["Pos"]))

    # ============================================================================
    # INTEGRATION TESTS - FULL PIPELINE
    # ============================================================================

    def test_ml_predictions_output_basic(self):
        """Test that ML predictions produce expected output files with correct structure.
        This test saves output to a fixed directory for inspection."""
        # Use a fixed output directory for inspection
        output_dir = os.path.join(self.test_data_path, "ml_predictor_test_output")
        os.makedirs(output_dir, exist_ok=True)
        
        mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
        mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
        mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
        
        result = run_ml_predictions(
            mhc1_agg_file,
            mhc1_all_file,
            mhc2_agg_file,
            self.model_artifacts_path,
            output_dir,
            'HCC1395',
            0.55
        )
        
        # Check that output file was created
        self.assertTrue(os.path.exists(result))
        
        # Check file path matches expected naming convention
        self.assertIn("_predict_pvacview.tsv", result)
        
        # Check that file contains data
        df = pd.read_csv(result, sep='\t')
        
        self.assertGreater(df.shape[0], 0)
        
        # Check required columns exist
        self.assertIn('Evaluation', df.columns)
        self.assertIn('ML Prediction (score)', df.columns)
        
        # Check Evaluation values are valid
        valid_evaluations = {'Accept', 'Reject', 'Pending'}
        self.assertTrue(df['Evaluation'].isin(valid_evaluations).all())
        
        # Print file path for easy inspection
        print(f"\nTest output file saved to: {result}")
        
        #df.to_csv(os.path.join(output_dir, "HCC1395_predict_pvacview.tsv"), sep='\t', index=False)

    def test_ml_predictions_different_threshold(self):
        """Test that different threshold values produce different results."""
        with tempfile.TemporaryDirectory() as output_dir_low:
            with tempfile.TemporaryDirectory() as output_dir_high:
                mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
                mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
                mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
                
                # Run with low threshold
                result_low = run_ml_predictions(
                    mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                    self.model_artifacts_path, output_dir_low, 'HCC1395', 0.30
                )
                
                # Run with high threshold
                result_high = run_ml_predictions(
                    mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                    self.model_artifacts_path, output_dir_high, 'HCC1395', 0.80
                )
                
                df_low = pd.read_csv(result_low, sep='\t')
                df_high = pd.read_csv(result_high, sep='\t')
                
                # With lower threshold, we should have more Accept predictions
                # (assuming the model produces probabilities across the range)
                self.assertGreaterEqual(
                    (df_low['Evaluation'] == 'Accept').sum(),
                    (df_high['Evaluation'] == 'Accept').sum()
                )


    def test_output_format_new_format_version(self):
        """Test the format of the output file."""
        with tempfile.TemporaryDirectory() as output_dir:
            mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
            mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
            mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
            
            result = run_ml_predictions(
                mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                self.model_artifacts_path, output_dir, 'HCC1395', 0.55
            )
            
            df = pd.read_csv(result, sep='\t')
            
            # Check ML Prediction (score) format
            accept_rows = df[df['Evaluation'] == 'Accept']
            if len(accept_rows) > 0:
                # Accept rows should have format "Accept (X.XX)"
                valid_scores = accept_rows['ML Prediction (score)'].dropna()
                if len(valid_scores) > 0:
                    self.assertTrue(
                        valid_scores.str.startswith("Accept (").all(),
                        "Some Accept rows don't have correct ML Prediction (score) format"
                    )
            
            pending_rows = df[df['Evaluation'] == 'Pending']
            if len(pending_rows) > 0:
                # Pending rows can have either "NA" or "Review (X.XX)" format
                # (because Review gets converted to Pending but keeps original score format)
                pending_scores = pending_rows['ML Prediction (score)'].dropna()
                if len(pending_scores) > 0:
                    expected_formats = (
                        pending_scores == "NA"
                    ) | (
                        pending_scores.str.startswith("Review (")
                    )
                    self.assertTrue(
                        expected_formats.all(),
                        "Some Pending rows don't have correct ML Prediction (score) format"
                    )

    def test_output_preserves_original_columns(self):
        """Test that output files preserve columns from original aggregated file."""
        with tempfile.TemporaryDirectory() as output_dir:
            mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
            mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
            mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
            
            # Read original file to get expected columns
            original_df = pd.read_csv(mhc1_agg_file, sep='\t', dtype=str)
            original_cols = set(original_df.columns)
            original_cols.discard('Evaluation')  # This will be replaced
            
            result = run_ml_predictions(
                mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                self.model_artifacts_path, output_dir, 'HCC1395', 0.55
            )
            
            output_df = pd.read_csv(result, sep='\t')
            output_cols = set(output_df.columns)
            output_cols.discard('Comments')  # This is new
            
            # All original columns (except Evaluation) should be present
            self.assertTrue(original_cols.issubset(output_cols))

    def test_tsl_column_format(self):
        """Test that TSL column is properly formatted as integer in output."""
        with tempfile.TemporaryDirectory() as output_dir:
            mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
            mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
            mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
            
            result = run_ml_predictions(
                mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                self.model_artifacts_path, output_dir, 'HCC1395', 0.55
            )
            
            df = pd.read_csv(result, sep='\t')
            
            # TSL should be integer if it exists
            if 'TSL' in df.columns:
                self.assertTrue(pd.api.types.is_integer_dtype(df['TSL'].fillna(6)))

    # ============================================================================
    # EDGE CASE TESTS
    # ============================================================================

    def test_merge_with_mismatched_row_counts(self):
        """Test behavior when MHC Class I and Class II files have different row counts."""
        # This tests the warning path in merge_and_prepare_data
        mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
        mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
        mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
        
        # The function should handle this gracefully
        result = merge_and_prepare_data(mhc1_agg_file, mhc1_all_file, mhc2_agg_file)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(result.shape[0], 0)

    def test_pending_evaluation_for_unmatched_ids(self):
        """Test that unmatched IDs get 'Pending' evaluation."""
        with tempfile.TemporaryDirectory() as output_dir:
            mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
            mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
            mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
            
            result = run_ml_predictions(
                mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                self.model_artifacts_path, output_dir, 'HCC1395', 0.55
            )
            
            df = pd.read_csv(result, sep='\t')
            
            # If there are any NA evaluations, they should be set to Pending
            self.assertFalse(df['Evaluation'].isna().any())

    # ============================================================================
    # ERROR HANDLING TESTS
    # ============================================================================

    def test_missing_input_file_error(self):
        """Test that missing input files raise appropriate errors."""
        with tempfile.TemporaryDirectory() as output_dir:
            mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "nonexistent_file.tsv")
            mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
            mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
            
            with self.assertRaises(Exception):
                run_ml_predictions(
                    mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                    self.model_artifacts_path, output_dir, 'HCC1395', 0.55
                )

    def test_missing_model_artifacts_error(self):
        """Test that missing model artifacts raise appropriate errors."""
        with tempfile.TemporaryDirectory() as output_dir:
            with tempfile.TemporaryDirectory() as fake_artifacts_dir:
                mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
                mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
                mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
                
                with self.assertRaises(Exception):
                    run_ml_predictions(
                        mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                        fake_artifacts_dir, output_dir, 'HCC1395', 0.55
                    )

    # ============================================================================
    # DATA VALIDATION TESTS
    # ============================================================================

    def test_probability_values_valid_range(self):
        """Test that prediction probabilities are in valid range [0, 1]."""
        with tempfile.TemporaryDirectory() as output_dir:
            mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
            mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
            mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
            
            result = run_ml_predictions(
                mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                self.model_artifacts_path, output_dir, 'HCC1395', 0.55
            )
            
            # Read the output file and extract probabilities
            df = pd.read_csv(result, sep='\t')
            accept_rows = df[df['Evaluation'] == 'Accept']
            
            if len(accept_rows) > 0:
                # Extract probability from ML Prediction (score) column
                scores = accept_rows['ML Prediction (score)'].str.extract(r'Accept \(([\d.]+)\)')
                probabilities = pd.to_numeric(scores[0])
                self.assertTrue((probabilities >= 0).all())
                self.assertTrue((probabilities <= 1).all())

    def test_output_row_count_matches_input(self):
        """Test that output has same number of rows as input aggregated file."""
        with tempfile.TemporaryDirectory() as output_dir:
            mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
            mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
            mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
            
            original_df = pd.read_csv(mhc1_agg_file, sep='\t')
            original_row_count = len(original_df)
            
            result = run_ml_predictions(
                mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                self.model_artifacts_path, output_dir, 'HCC1395', 0.55
            )
            
            output_df = pd.read_csv(result, sep='\t')
            output_row_count = len(output_df)
            
            # Output should have same number of rows as input
            self.assertEqual(original_row_count, output_row_count)

    # ============================================================================
    # SAMPLE NAME TESTS
    # ============================================================================

    def test_output_file_naming(self):
        """Test that output files use correct sample name."""
        with tempfile.TemporaryDirectory() as output_dir:
            mhc1_agg_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv")
            mhc1_all_file = os.path.join(self.test_data_path, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.tsv")
            mhc2_agg_file = os.path.join(self.test_data_path, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv")
            
            sample_name = "TEST_SAMPLE"
            result = run_ml_predictions(
                mhc1_agg_file, mhc1_all_file, mhc2_agg_file,
                self.model_artifacts_path, output_dir, sample_name, 0.55
            )
            
            self.assertIn(sample_name, result)


if __name__ == '__main__':
    unittest.main()

