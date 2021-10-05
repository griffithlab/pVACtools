import unittest
import unittest.mock
import os
import sys
import re
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile
import pandas as pd
from mock import patch

import pvactools.lib.call_iedb
from pvactools.lib.prediction_class import PredictionClass, IEDB

from .test_utils import *

def test_data_directory():
    return os.path.join(pvactools_directory(), 'tests', 'test_data', 'call_iedb')

class CallIEDBTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'call_iedb.py')
        cls.test_data_dir  = test_data_directory()
        cls.additional_setup()

    @classmethod
    def additional_setup(cls):
        pass

class CallIEDBCompileTests(CallIEDBTests):
    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

class FilterResponseTests(CallIEDBTests):
    def test_filter_response_ok(self):
        unfiltered_file = os.path.join(self.test_data_dir, 'unfiltered.txt')
        filtered_file   = os.path.join(self.test_data_dir, 'filtered.txt')
        with open(unfiltered_file, 'rb') as f:
            unfiltered_file_contents = f.read().rstrip()
        with open(filtered_file, 'rb') as f:
            filtered_file_contents = f.read().rstrip()
        filtered_response = IEDB.filter_response(unfiltered_file_contents)
        self.assertEqual(filtered_response, filtered_file_contents)
        filtered_response_on_filtered_file = IEDB.filter_response(filtered_file_contents)
        self.assertEqual(filtered_response_on_filtered_file, filtered_file_contents)

class CallIEDBClassITests(CallIEDBTests):
    @classmethod
    def additional_setup(cls):
        cls.input_file     = os.path.join(cls.test_data_dir, 'input.fasta')
        cls.allele         = 'HLA-A*02:01'
        cls.epitope_length = 9
        cls.methods = ['ann', 'smmpmbec', 'smm']

    def test_iedb_methods_generate_expected_files(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request:
            #netmhcpan, netmhccons, and pickpocket are slow so we won't run them in the tests
            for method in self.methods:
                call_iedb_output_file = tempfile.NamedTemporaryFile()
                class_name = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)

                pvactools.lib.call_iedb.main([
                    self.input_file,
                    call_iedb_output_file.name,
                    class_name,
                    self.allele,
                    '-l', str(self.epitope_length)
                ])
                mock_request.assert_has_calls([
                    generate_class_i_call(method, self.allele, self.epitope_length, self.input_file)
                ])
                reader = open(self.input_file, mode='r')
                reader.close()
                expected_output_file = os.path.join(self.test_data_dir, 'output_%s.tsv' % method)
                self.assertTrue(cmp(call_iedb_output_file.name, expected_output_file))

    #the output from MHCflurry varies between operating systems and the version of tensorflow installed
    #these outputs where created on tensorflow 1.8.0
    def test_mhcflurry_method_generates_expected_files(self):
        call_iedb_output_file = tempfile.NamedTemporaryFile()

        pvactools.lib.call_iedb.main([
            self.input_file,
            call_iedb_output_file.name,
            'MHCflurry',
            self.allele,
            '-l', str(self.epitope_length)
        ])
        if sys.platform == 'darwin':
            expected_output_file = os.path.join(self.test_data_dir, 'output_mhcflurry_osx.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[1,6,7])
            actual_df = pd.read_csv(call_iedb_output_file.name, sep="\t", index_col=[1,6,7])
            pd.testing.assert_frame_equal(expected_df, actual_df, check_like=True, check_less_precise=0)

    def test_mhcnuggets_method_generates_expected_files(self):
        call_iedb_output_file = tempfile.NamedTemporaryFile()

        pvactools.lib.call_iedb.main([
            self.input_file,
            call_iedb_output_file.name,
            'MHCnuggetsI',
            self.allele,
            '-l', str(self.epitope_length)
        ])
        expected_output_file = os.path.join(self.test_data_dir, 'output_mhcnuggetsI.tsv')
        expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,2,3])
        actual_df = pd.read_csv(call_iedb_output_file.name, sep="\t", index_col=[0,2,3])
        pd.testing.assert_frame_equal(expected_df, actual_df, check_like=True, check_less_precise=0)

class CallIEDBClassIITests(CallIEDBTests):
    @classmethod
    def additional_setup(cls):
        cls.input_file     = os.path.join(cls.test_data_dir, 'input_31.fasta')
        cls.allele         = 'H2-IAb'
        cls.methods = ['nn_align']
        cls.request_mock = unittest.mock.Mock(side_effect = (
            make_response(method, cls.test_data_dir) for method in cls.methods
        ))
        pvactools.lib.call_iedb.requests.post = cls.request_mock

    def test_iedb_methods_generate_expected_files(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request:
            for method in self.methods:
                call_iedb_output_file = tempfile.NamedTemporaryFile()
                class_name = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)

                pvactools.lib.call_iedb.main([
                    self.input_file,
                    call_iedb_output_file.name,
                    class_name,
                    self.allele,
                    '-l', '15',
                ])
                mock_request.assert_has_calls([
                    generate_class_ii_call(method, self.allele, 15, self.input_file)
                ])
                expected_output_file = os.path.join(self.test_data_dir, 'output_%s.tsv' % method)
                self.assertTrue(cmp(call_iedb_output_file.name, expected_output_file))

    def test_mhcnuggets_method_generates_expected_files(self):
        call_iedb_output_file = tempfile.NamedTemporaryFile()

        pvactools.lib.call_iedb.main([
            self.input_file,
            call_iedb_output_file.name,
            'MHCnuggetsII',
            'DPA1*01:03',
            '-l', '15',
        ])
        expected_output_file = os.path.join(self.test_data_dir, 'output_mhcnuggetsII.tsv')
        expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,2,3])
        actual_df = pd.read_csv(call_iedb_output_file.name, sep="\t", index_col=[0,2,3])
        pd.testing.assert_frame_equal(expected_df, actual_df, check_like=True, check_less_precise=0)

if __name__ == '__main__':
    unittest.main()
