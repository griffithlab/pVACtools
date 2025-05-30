import unittest
import os
import tempfile
import json
from pvactools.tools.pvaccompare.runners.run_compare_aggregated_tsv import main


# python -m unittest tests/test_compare_aggregated_tsv.py
# python -m unittest discover -s tests
class TestRunCompareAggregatedTSV(unittest.TestCase):
    def setUp(self):
        self.input_file1 = tempfile.NamedTemporaryFile(delete=False, suffix=".json")
        self.input_file2 = tempfile.NamedTemporaryFile(delete=False, suffix=".json")
        self.temp_dir = tempfile.TemporaryDirectory()
        self.columns_to_compare = [
            "Num Passing Transcripts",
            "Best Peptide",
            "Best Transcript",
            "Num Passing Peptides",
            "Tier",
        ]
        self.output_path = self.temp_dir.name
        self.file_name = "aggregated_data.json"
        self.data_dir = "tests/test_data/pvaccompare/"
        self.class_type = "2"

    def tearDown(self):
        os.remove(self.input_file1.name)
        os.remove(self.input_file2.name)
        self.temp_dir.cleanup()

    def test_identical_files(self):
        with open(self.data_dir + "aggregated_input1.tsv", "r") as f:
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

        self.assertIn("INFO:root:The Aggregated TSV files are identical.", log.output)

    def test_different_files(self):
        expected_output_path = self.data_dir + "aggregated_expected_output.json"

        with open(self.data_dir + "aggregated_input1.tsv", "r") as f:
            content1 = f.read()
        with open(self.data_dir + "aggregated_input2.tsv", "r") as f:
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
            "INFO:root:• Renamed 'best peptide' to 'Best Peptide' in file 1", log.output
        )
        self.assertIn(
            "INFO:root:• Renamed 'best transcript' to 'Best Transcript' in file 2",
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

    def test_missing_id(self):
        expected_output_path = self.data_dir + "aggregated_id_change_output.json"

        with open(self.data_dir + "aggregated_input1.tsv", "r") as f:
            content1 = f.read()
        with open(self.data_dir + "aggregated_input3.tsv", "r") as f:
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
        self.assertIn("INFO:root:• Replaced ID with Gene and AA Change", log.output)

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
