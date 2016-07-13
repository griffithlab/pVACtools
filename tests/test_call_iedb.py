import unittest
import unittest.mock
import os
import sys
import re
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile
base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
sys.path.append(base_dir)
import pvacseq.lib.call_iedb

def unordered_test(file1, file2):
    r1 = open(file1, mode='r')
    s1 = set(row.rstrip() for row in r1)
    r1.close()
    r2 = open(file2, mode='r')
    s2 = set(row.rstrip() for row in r2)
    r2.close()
    result = s1^s2
    print(result)
    return not len(result)

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
        cls.executable_dir = os.path.join(base_dir, 'pvacseq', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'call_iedb.py')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'call_iedb')
        cls.input_file     = os.path.join(cls.test_data_dir, 'input.fasta')
        cls.allele         = 'HLA-A02:01'
        cls.epitope_length = 9
        cls.methods = ['ann', 'smmpmbec', 'smm']
        cls.request_mock = unittest.mock.Mock(side_effect = (
            make_response(method, cls.test_data_dir) for method in cls.methods
        ))
        pvacseq.lib.call_iedb.requests.post = cls.request_mock


    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_iedb_methods_generate_expected_files(self):
        #netmhcpan, netmhccons, and pickpocket are slow so we won't run them in the tests
        for method in self.methods:
            # call_iedb_output_file = tempfile.NamedTemporaryFile()
            call_iedb_output_file = lambda :None
            call_iedb_output_file.name = "testout_%s.tsv"%method

            call_iedb_command = "%s %s %s %s %s %s %s" % (
                self.python,
                self.executable,
                self.input_file,
                call_iedb_output_file.name,
                method,
                self.allele,
                self.epitope_length,
            )

            # self.assertFalse(call(call_iedb_command, shell=True))
            pvacseq.lib.call_iedb.main([
                self.input_file,
                call_iedb_output_file.name,
                method,
                self.allele,
                str(self.epitope_length)
            ])
            reader = open(self.input_file, mode='r')
            self.request_mock.assert_called_with('http://tools-api.iedb.org/tools_api/mhci/', data={
                'sequence_text':reader.read(),
                'method': method,
                'allele': re.sub(r'(\w*-[\w|\d])(.*)', r'\1*\2', self.allele),
                'length': self.epitope_length
            })
            reader.close()
            expected_output_file = os.path.join(self.test_data_dir, 'output_%s.tsv' % method)
            self.assertTrue(cmp(call_iedb_output_file.name, expected_output_file))
            # self.assertTrue(unordered_test(call_iedb_output_file.name, expected_output_file))

if __name__ == '__main__':
    unittest.main()
