import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

class CallIEDBTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'pvacseq', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'call_iedb.py')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'call_iedb')
        cls.input_file     = os.path.join(cls.test_data_dir, 'input.fasta')
        cls.allele         = 'HLA-A02:01'
        cls.epitope_length = 9

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_iedb_methods_generate_expected_files(self):
        #netmhcpan, netmhccons, and pickpocket are slow so we won't run them in the tests
        methods = ['ann', 'smmpmbec', 'smm']
        for method in methods:
            call_iedb_output_file = tempfile.NamedTemporaryFile()

            call_iedb_command = "%s %s %s %s %s %s %s" % (
                self.python,
                self.executable,
                self.input_file,
                call_iedb_output_file.name,
                method,
                self.allele,
                self.epitope_length,
            )

            self.assertFalse(call(call_iedb_command, shell=True))
            expected_output_file = os.path.join(self.test_data_dir, 'output_%s.tsv' % method)
            self.assertTrue(cmp(call_iedb_output_file.name, expected_output_file))

if __name__ == '__main__':
    unittest.main()
