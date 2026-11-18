import unittest
import os
import requests
import pandas as pd
import json
import time

from tests.utils import *

class IEDBApiStatusTest(unittest.TestCase):
    def test_iedb_class_i_api_runs_and_returns_expected_value(self):
        data = {
          "pipeline_id": "",
          "run_stage_range": [1,1],
          "stages": [
            {
              "stage_number": 1,
              "tool_group": "mhci",
              "input_sequence_text": ">1\nASTPGHTIIYEAVCLHNDRTTIP",
              "input_parameters":
                {
                  "alleles": "HLA-A*02:01",
                  "peptide_length_range": [9,9],
                  "predictors": [
                    {
                    "type": "binding",
                    "method": "ann",
                    }
                  ]
                }
            }
          ]
        }

        query_response = requests.post("https://api-nextgen-tools.iedb.org/api/v1/pipeline", json=data)
        self.assertEqual(query_response.status_code, 200)

        results_uri = json.loads(query_response.text)['results_uri']
        results_response = requests.get(results_uri)
        self.assertEqual(results_response.status_code, 200)
        results = json.loads(results_response.text)
        while (results['status'] == 'pending'):
            time.sleep(5)
            results_response = requests.get(results_uri)
            self.assertEqual(results_response.status_code, 200)
            results = json.loads(results_response.text)

        self.assertEqual(len(results['data']['errors']), 0)
        result_table = [i for i in results['data']['results'] if i['type'] == 'peptide_table'][0]
        data = result_table['table_data']
        headers = [c['name'].replace(f"ann_","") for c in result_table['table_columns']]
        results_df = pd.DataFrame(data, columns=headers)
        results_df.rename(columns={
            'peptide_index': 'seq_num',
        }, inplace=True)
        expected_file_name = os.path.join(pvactools_directory(), 'api_status_tests', 'iedb_class_i_response.tsv')
        expected_df = pd.read_csv(expected_file_name, sep="\t")
        pd.testing.assert_frame_equal(results_df, expected_df, check_dtype=False)

    def test_iedb_class_ii_api_runs_and_returns_expected_value(self):
        data = {
          "pipeline_id": "",
          "run_stage_range": [1,1],
          "stages": [
            {
              "stage_number": 1,
              "tool_group": "mhcii",
              "input_sequence_text": ">1\nASTPGHTIIYEAVCLHNDRTTIP",
              "input_parameters":
                {
                  "alleles": "HLA-DPA1*01:03/DPB1*02:01",
                  "peptide_length_range": [15,15],
                  "predictors": [
                    {
                    "type": "binding",
                    "method": "nn_align",
                    }
                  ]
                }
            }
          ]
        }

        query_response = requests.post("https://api-nextgen-tools.iedb.org/api/v1/pipeline", json=data)
        self.assertEqual(query_response.status_code, 200)

        results_uri = json.loads(query_response.text)['results_uri']
        results_response = requests.get(results_uri)
        self.assertEqual(results_response.status_code, 200)
        results = json.loads(results_response.text)
        while (results['status'] == 'pending'):
            time.sleep(5)
            results_response = requests.get(results_uri)
            self.assertEqual(results_response.status_code, 200)
            results = json.loads(results_response.text)

        self.assertEqual(len(results['data']['errors']), 0)
        result_table = [i for i in results['data']['results'] if i['type'] == 'peptide_table'][0]
        data = result_table['table_data']
        headers = [c['name'].replace(f"ann_","") for c in result_table['table_columns']]
        results_df = pd.DataFrame(data, columns=headers)
        results_df.rename(columns={
            'peptide_index': 'seq_num',
        }, inplace=True)
        expected_file_name = os.path.join(pvactools_directory(), 'api_status_tests', 'iedb_class_ii_response.tsv')
        expected_df = pd.read_csv(expected_file_name, sep="\t")
        pd.testing.assert_frame_equal(results_df, expected_df, check_dtype=False)
