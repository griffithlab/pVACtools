import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

class ParseOutputTests(unittest.TestCase):
    def setUp(self):
        self.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        self.executable_dir = os.path.join(base_dir, 'pvacseq')
        self.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'parse_output')
        self.sample_name             = 'Test'
        self.peptide_sequence_length = 21
        self.epitope_length          = 9
        self.allele                  = 'HLA-A29:02'

    def tearDown(self):
        del self.python
        del self.executable_dir
        del self.test_data_dir
        del self.sample_name
        del self.peptide_sequence_length
        del self.epitope_length
        del self.allele

    def test_parse_output_runs_and_produces_expected_output(self):
        parse_output_netmhc_executable  = os.path.join(self.executable_dir, 'parse_output.py')
        self.assertTrue(py_compile.compile(parse_output_netmhc_executable))

        parse_output_netmhc_input_file  = os.path.join(self.test_data_dir, "%s.%s.%s.netmhc.xls" % (self.sample_name, self.allele, self.epitope_length))
        parse_output_netmhc_key_file  = os.path.join(self.test_data_dir, "%s_%s.key" % (self.sample_name, self.peptide_sequence_length))
        parse_output_netmhc_output_file = tempfile.NamedTemporaryFile().name

        parse_output_netmhc_command = "%s %s %s %s %s" % (self.python, parse_output_netmhc_executable, parse_output_netmhc_input_file, parse_output_netmhc_key_file, parse_output_netmhc_output_file)
        self.assertFalse(call(parse_output_netmhc_command, shell=True))

        expected_output_file  = os.path.join(self.test_data_dir, "%s.%s.%s.netmhc.parsed" % (self.sample_name, self.allele, self.epitope_length))
        self.assertTrue(cmp(parse_output_netmhc_output_file, expected_output_file))
