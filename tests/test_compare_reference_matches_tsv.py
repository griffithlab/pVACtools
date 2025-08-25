import unittest
import os
import tempfile
import json
from pvactools.tools.pvaccompare.runners.run_compare_reference_matches_tsv import main


# python -m unittest tests/test_compare_reference_matches_tsv.py
# python -m unittest discover -s tests
class TestRunCompareReferenceMatchesTSV(unittest.TestCase):
    def setUp(self):
        self.input_file1 = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
        self.input_file2 = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
        self.temp_dir = tempfile.TemporaryDirectory()
        self.columns_to_compare = ["Peptide", "Match Window"]
        self.output_path = self.temp_dir.name
        self.file_name = "reference_matches_data.json"
        self.data_dir = "tests/test_data/pvaccompare/"
        self.class_type = "1"

    def tearDown(self):
        os.remove(self.input_file1.name)
        os.remove(self.input_file2.name)
        self.temp_dir.cleanup()

    def test_identical_files(self):
        with open(self.data_dir + "reference_matches_input1.tsv", "r") as f:
            content = f.read()
        self.input_file1.write(content.encode())
        self.input_file2.write(content.encode())
        self.input_file1.close()
        self.input_file2.close()

        with self.assertLogs(level="INFO") as log:
            main(
                self.input_file1.name,
                self.input_file2.name,
                self.columns_to_compare,
                self.output_path,
                self.class_type,
            )

        self.assertIn(
            "INFO:root:The Reference Matches TSV files are identical.", log.output
        )

    def test_different_files(self):
        expected_output_path = self.data_dir + "reference_matches_expected_output_normal.json"

        with open(self.data_dir + "reference_matches_input1.tsv", "r") as f:
            content1 = f.read()
        with open(self.data_dir + "reference_matches_input2.tsv", "r") as f:
            content2 = f.read()

        self.input_file1.write(content1.encode())
        self.input_file2.write(content2.encode())
        self.input_file1.close()
        self.input_file2.close()

        main(
            self.input_file1.name,
            self.input_file2.name,
            self.columns_to_compare,
            self.output_path,
            self.class_type,
        )

        with open(f"{self.output_path}/{self.file_name}") as f1, open(
            expected_output_path
        ) as f2:
            output_json_data = json.load(f1)
            expected_output = json.load(f2)

        output_json_data.pop("input_file1", None)
        output_json_data.pop("input_file2", None)
        expected_output.pop("input_file1", None)
        expected_output.pop("input_file2", None)

        self.assertEqual(output_json_data, expected_output)

    def test_duplicate_records(self):
        expected_output_path = self.data_dir + "reference_matches_expected_output_duplicates.json"

        with open(self.data_dir + "reference_matches_input1.tsv", "r") as f:
            content1 = f.read()
        with open(self.data_dir + "reference_matches_input3.tsv", "r") as f:
            content2 = f.read()

        self.input_file1.write(content1.encode())
        self.input_file2.write(content2.encode())
        self.input_file1.close()
        self.input_file2.close()

        with self.assertLogs(level="INFO") as log:
            main(
                self.input_file1.name,
                self.input_file2.name,
                self.columns_to_compare,
                self.output_path,
                self.class_type,
            )
        self.assertIn(
            "ERROR:root:ERROR: Duplicate unique records were found in file 2. Writing number of hits only.",
            log.output,
        )

        with open(f"{self.output_path}/{self.file_name}") as f1, open(
            expected_output_path
        ) as f2:
            output_json_data = json.load(f1)
            expected_output = json.load(f2)

        output_json_data.pop("input_file1", None)
        output_json_data.pop("input_file2", None)
        expected_output.pop("input_file1", None)
        expected_output.pop("input_file2", None)

        self.assertEqual(output_json_data, expected_output)
