import unittest
import os
import tempfile
import json
from pvactools.tools.pvaccompare.runners.run_compare_json import main


# python -m unittest tests/test_compare_json.py
# python -m unittest discover -s tests
class TestRunCompareJSON(unittest.TestCase):
    def setUp(self):
        self.input_file1 = tempfile.NamedTemporaryFile(delete=False, suffix=".json")
        self.input_file2 = tempfile.NamedTemporaryFile(delete=False, suffix=".json")
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_path = self.temp_dir.name
        self.file_name = "json_input_data.json"
        self.data_dir = "tests/test_data/pvaccompare/"
        self.class_type = "2"

    def tearDown(self):
        os.remove(self.input_file1.name)
        os.remove(self.input_file2.name)
        self.temp_dir.cleanup()

    def test_identical_files(self):
        with open(self.data_dir + "json_input1.json", "r") as f:
            content = f.read()
        self.input_file1.write(content.encode())
        self.input_file2.write(content.encode())
        self.input_file1.close()
        self.input_file2.close()

        with self.assertLogs(level="INFO") as log:
            main(
                self.input_file1.name,
                self.input_file2.name,
                self.output_path,
                self.class_type,
            )

        self.assertIn("INFO:root:The JSON metric inputs are identical.", log.output)

    def test_different_files(self):
        expected_output_path = self.data_dir + "json_expected_output.json"

        with open(self.data_dir + "json_input1.json", "r") as f:
            content1 = f.read()
        with open(self.data_dir + "json_input2.json", "r") as f:
            content2 = f.read()

        self.input_file1.write(content1.encode())
        self.input_file2.write(content2.encode())
        self.input_file1.close()
        self.input_file2.close()

        main(
            self.input_file1.name,
            self.input_file2.name,
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
    
    def test_different_alleles(self):
        expected_output_path = self.data_dir + "json_diff_alleles_expected_output.json"

        with open(self.data_dir + "json_input1.json", "r") as f:
            content1 = f.read()
        with open(self.data_dir + "json_input3.json", "r") as f:
            content2 = f.read()

        self.input_file1.write(content1.encode())
        self.input_file2.write(content2.encode())
        self.input_file1.close()
        self.input_file2.close()

        main(
            self.input_file1.name,
            self.input_file2.name,
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
