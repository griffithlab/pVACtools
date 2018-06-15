import unittest
import unittest.mock
import os
import sys
import re
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile
import lib.call_iedb
from lib.prediction_class import PredictionClass, IEDB

def make_response(method, path):
    reader = open(os.path.join(
        path,
        'response_%s.tsv' %method
    ), mode='r')
    response_obj = lambda :None
    response_obj.status_code = 200
    response_obj.text = reader.read()
    reader.close()
    return response_obj

class CallIEDBTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'call_iedb.py')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'call_iedb')
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
        cls.short_input_file = os.path.join(cls.test_data_dir, 'input.short.fasta')
        cls.allele         = 'HLA-A*02:01'
        cls.epitope_length = 9
        cls.methods = ['ann', 'smmpmbec', 'smm']
        cls.request_mock = unittest.mock.Mock(side_effect = (
            make_response(method, cls.test_data_dir) for method in cls.methods
        ))
        lib.call_iedb.requests.post = cls.request_mock

    def test_iedb_methods_generate_expected_files(self):
        #netmhcpan, netmhccons, and pickpocket are slow so we won't run them in the tests
        for method in self.methods:
            call_iedb_output_file = tempfile.NamedTemporaryFile()
            class_name = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)

            lib.call_iedb.main([
                self.input_file,
                call_iedb_output_file.name,
                class_name,
                self.allele,
                '-l', str(self.epitope_length)
            ])
            reader = open(self.input_file, mode='r')
            self.request_mock.assert_called_with('http://tools-cluster-interface.iedb.org/tools_api/mhci/', data={
                'sequence_text':reader.read(),
                'method': method,
                'allele': self.allele,
                'length': self.epitope_length,
                'user_tool': 'pVac-seq',
            })
            reader.close()
            expected_output_file = os.path.join(self.test_data_dir, 'output_%s.tsv' % method)
            self.assertTrue(cmp(call_iedb_output_file.name, expected_output_file))

    #the output from MHCflurry varies between operating systems and the version of tensorflow installed
    #these outputs where created on tensorflow 1.1.0
    def test_mhcflurry_method_generates_expected_files(self):
        call_iedb_output_file = tempfile.NamedTemporaryFile()

        lib.call_iedb.main([
            self.input_file,
            call_iedb_output_file.name,
            'MHCflurry',
            self.allele,
            '-l', str(self.epitope_length)
        ])
        if sys.platform == 'darwin':
            expected_output_file = os.path.join(self.test_data_dir, 'output_mhcflurry_osx.tsv')
            self.assertTrue(cmp(call_iedb_output_file.name, expected_output_file))

    def test_mhcnuggets_method_generates_expected_files(self):
        call_iedb_output_file = tempfile.NamedTemporaryFile()

        lib.call_iedb.main([
            self.short_input_file,
            call_iedb_output_file.name,
            'MHCnuggetsI',
            'HLA-A0201',
            '-l', str(self.epitope_length)
        ])
        expected_output_file = os.path.join(self.test_data_dir, 'output_mhcnuggetsI.tsv')
        self.assertTrue(cmp(call_iedb_output_file.name, expected_output_file))

class CallIEDBClassIITests(CallIEDBTests):
    @classmethod
    def additional_setup(cls):
        cls.input_file     = os.path.join(cls.test_data_dir, 'input_31.fasta')
        cls.short_input_file = os.path.join(cls.test_data_dir, 'input_31.short.fasta')
        cls.allele         = 'H2-IAb'
        cls.methods = ['nn_align']
        cls.request_mock = unittest.mock.Mock(side_effect = (
            make_response(method, cls.test_data_dir) for method in cls.methods
        ))
        lib.call_iedb.requests.post = cls.request_mock

    def test_iedb_methods_generate_expected_files(self):
        for method in self.methods:
            call_iedb_output_file = tempfile.NamedTemporaryFile()
            class_name = PredictionClass.prediction_class_name_for_iedb_prediction_method(method)

            lib.call_iedb.main([
                self.input_file,
                call_iedb_output_file.name,
                class_name,
                self.allele,
            ])
            reader = open(self.input_file, mode='r')
            self.request_mock.assert_called_with('http://tools-cluster-interface.iedb.org/tools_api/mhcii/', data={
                'sequence_text':reader.read(),
                'method': method,
                'allele': self.allele,
                'user_tool': 'pVac-seq',
            })
            reader.close()
            expected_output_file = os.path.join(self.test_data_dir, 'output_%s.tsv' % method)
            self.assertTrue(cmp(call_iedb_output_file.name, expected_output_file))

    def test_mhcnuggets_method_generates_expected_files(self):
        call_iedb_output_file = tempfile.NamedTemporaryFile()

        lib.call_iedb.main([
            self.short_input_file,
            call_iedb_output_file.name,
            'MHCnuggetsII',
            'H-2-IAb',
        ])
        expected_output_file = os.path.join(self.test_data_dir, 'output_mhcnuggetsII.tsv')
        self.assertTrue(cmp(call_iedb_output_file.name, expected_output_file))

if __name__ == '__main__':
    unittest.main()
