import unittest
import subprocess
import requests
import time
import py_compile
import argparse
import os
import tempfile
import shutil
import argparse
from unittest.mock import MagicMock
import pvactools.tools.compare as compare
import pvactools.tools.pvaccompare.compare_tools.validators as comparison_validators

# python -m unittest tests/test_compare.py
# python -m unittest discover -s tests
class TestRunCompare(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.url = "http://localhost:8080/main.html"
        cls.data_dir = "tests/test_data/pvaccompare/"
        cls.aggregated_columns = ["Best Peptide", "Best Transcript"]
        cls.unaggregated_columns = ["Median MT IC50 Score", "Median WT IC50 Score"]
        cls.reference_match_columns = ["Peptide", "Hit Definition"]
        cls.server_process = subprocess.Popen(
            ["python", "-m", "pvactools.tools.pvaccompare.server"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        time.sleep(1)


    @classmethod
    def tearDownClass(cls):
        cls.server_process.terminate()
        cls.server_process.wait()


    def test_compare_server(self):
        response = requests.get(self.url)
        self.assertEqual(response.status_code, 200)


    def test_cors_headers(self):
        response = requests.get(self.url)
        self.assertEqual(response.headers.get("Access-Control-Allow-Origin"), "*")
        self.assertEqual(response.headers.get("Access-Control-Allow-Methods"), "GET, OPTIONS")
        self.assertIn("x-requested-with", response.headers.get("Access-Control-Allow-Headers", ""))


    def test_options_request(self):
        response = requests.options(self.url)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.headers.get("Access-Control-Allow-Methods"), "GET, OPTIONS")


    def test_compare_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            'pvactools',
            'tools',
            "compare.py"
        ))
        self.assertTrue(compiled_pvac_path)


    def test_compare_parser(self):
        parser = compare.define_parser()
        self.assertEqual(type(parser), argparse.ArgumentParser)


    def test_compare_aggregated_validator(self):
        aggregated_columns = self.aggregated_columns + ["Test1"]

        parser = compare.define_parser()
        parser.error = MagicMock()
        comparison_validators.validate_aggregated_columns(aggregated_columns, parser)
        parser.error.assert_called_once()
        error_msg = parser.error.call_args[0][0]
        self.assertIn("Invalid aggregated column 'Test1' specified.", error_msg)
    

    def test_compare_unaggregated_validator(self):
        unaggregated_columns = self.unaggregated_columns + ["Test1"]

        parser = compare.define_parser()
        parser.error = MagicMock()
        comparison_validators.validate_unaggregated_columns(unaggregated_columns, parser)
        parser.error.assert_called_once()
        error_msg = parser.error.call_args[0][0]
        self.assertIn("Invalid unaggregated column 'Test1' specified.", error_msg)
    

    def test_compare_reference_match_validator(self):
        reference_match_columns = self.reference_match_columns + ["Test1"]

        parser = compare.define_parser()
        parser.error = MagicMock()
        comparison_validators.validate_reference_match_columns(reference_match_columns, parser)
        parser.error.assert_called_once()
        error_msg = parser.error.call_args[0][0]
        self.assertIn("Invalid reference match column 'Test1' specified.", error_msg)


    def test_run_comparison(self):
        with tempfile.TemporaryDirectory() as input1, tempfile.TemporaryDirectory() as input2, tempfile.TemporaryDirectory() as output:
            input1_mhc_class_i = os.path.join(input1, "MHC_Class_I")
            input2_mhc_class_i = os.path.join(input2, "MHC_Class_I")
            output_mhc_class_i = os.path.join(output, "mhc_class_i")
            
            os.makedirs(input1_mhc_class_i)
            os.makedirs(input2_mhc_class_i)
            os.makedirs(output_mhc_class_i)

            test_data_files = {
                "all_epitopes.aggregated.metrics.json": self.data_dir + "json_input1.json",
                "all_epitopes.aggregated.tsv": self.data_dir + "aggregated_input1.tsv",
                "all_epitopes.tsv": self.data_dir + "unaggregated_input1.tsv",
                "reference_matches": self.data_dir + "reference_matches_input1.tsv"
            }

            for file_name, test_data_path in test_data_files.items():
                for folder in [input1_mhc_class_i, input2_mhc_class_i]:
                    file_path = os.path.join(folder, file_name)
                    shutil.copy(test_data_path, file_path)
            
            extra_all_epitopes_file = os.path.join(input2_mhc_class_i, "additional_all_epitopes.tsv")
            shutil.copy(test_data_files["all_epitopes.tsv"], extra_all_epitopes_file)

            args_list = [
                input1,
                input2,
                "--output-dir", output,
                "--aggregated-columns", "Best Peptide,Best Transcript",
                "--unaggregated-columns", "Median MT IC50 Score,Median WT IC50 Score",
                "--reference-match-columns", "Peptide,Hit Definition",
                "--mhc-class", "1",
                "--no-server"
            ]

            with self.assertLogs(level="INFO") as log:
                compare.main(args_list)

            self.assertIn("ERROR:root:ERROR: Could not locate the input YML file in results folder 1 for MHC Class I.", log.output)
            self.assertIn("ERROR:root:ERROR: Could not locate the input YML file in results folder 2 for MHC Class I.", log.output)
            self.assertIn("INFO:root:The JSON metric inputs are identical.", log.output)
            self.assertIn("ERROR:root:ERROR: Located multiple unaggregated TSV files in results folder 2 for MHC Class I:", log.output)
            self.assertIn("INFO:root:The Aggregated TSV files are identical.", log.output)
            self.assertIn("INFO:root:Successfully generated MHC Class I comparison report.", log.output)
