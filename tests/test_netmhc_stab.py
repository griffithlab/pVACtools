import unittest
import unittest.mock
from mock import patch
import logging
from testfixtures import LogCapture, StringComparison as S
import os
import re
import sys
import tempfile
import py_compile
from filecmp import cmp
import lib

def make_response(data, files, path, test_file):
    reader = open(os.path.join(
        path,
        test_file
    ), mode='rb')
    response_obj = lambda :None
    response_obj.status_code = 200
    response_obj.content = reader.read()
    reader.close()
    return response_obj

def make_rejected_response(data, files, path, self):
    if self.rejected:
        reader = open(os.path.join(
            path,
            'Netmhcstab.rejected.html'
        ), mode='rb')
        self.rejected = False
    else:
        reader = open(os.path.join(
            path,
            'Netmhcstab.html'
        ), mode='rb')
    response_obj = lambda :None
    response_obj.status_code = 200
    response_obj.content = reader.read()
    reader.close()
    return response_obj

class NetChopTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pvac_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.script_path = os.path.join(
            pvac_dir,
            'lib',
            'netmhc_stab.py'
        )
        cls.test_data_directory = os.path.join(
            pvac_dir,
            'tests',
            'test_data',
            'netmhc_stab'
        )
        cls.rejected = True

    def test_netmhc_stab_compiles(self):
        compiled_script_path = py_compile.compile(self.script_path)
        self.assertTrue(compiled_script_path)

    def test_netmhc_stab_runs(self):
        with patch('lib.netmhc_stab.requests.post',  unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            self.test_data_directory,
            'Netmhcstab.html'
           ))):
            output_file = tempfile.NamedTemporaryFile()
            lib.netmhc_stab.main([
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name
            ])
            self.assertTrue(cmp(
                os.path.join(self.test_data_directory, 'Test_filtered.stab.tsv'),
                output_file.name
            ))

    def test_netmhc_stab_fail(self):
        with patch('lib.netmhc_stab.requests.post',  unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            self.test_data_directory,
            'Netmhcstab.fail.html'
           ))), self.assertRaises(Exception) as context:
            output_file = tempfile.NamedTemporaryFile()
            lib.netmhc_stab.main([
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name
            ])
        self.assertTrue('NetMHCstabpan encountered an error during processing.' in str(context.exception))

    def test_netmhc_stab_rejected(self):
        logging.disable(logging.NOTSET)
        with patch('lib.netmhc_stab.requests.post',  unittest.mock.Mock(side_effect = lambda url, data, files=None: make_rejected_response(
            data,
           files,
           self.test_data_directory,
           self
           ))), LogCapture() as l:
            output_file = tempfile.NamedTemporaryFile()
            lib.netmhc_stab.main([
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name
            ])
            self.assertTrue(cmp(
                os.path.join(self.test_data_directory, 'Test_filtered.stab.tsv'),
                output_file.name
            ))
            l.check_present(('root', 'WARNING', S("Too many jobs submitted to NetMHCstabpan server. Waiting to retry.")))

    def test_netmhc_stab_cannot_open_file_error_retry(self):
        with patch('lib.netmhc_stab.requests.post',  unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            self.test_data_directory,
            'Netmhcstab.cannot_open_file_error.html'
           ))), self.assertRaises(Exception) as context:
            output_file = tempfile.NamedTemporaryFile()
            lib.netmhc_stab.main([
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name
            ])
        self.assertTrue('NetMHCstabpan server was unable to read the submitted fasta file.' in str(context.exception))

    #This is to ensure that we catch error cases that are not explicitly handled
    def test_netmhc_stab_other_error(self):
        with patch('lib.netmhc_stab.requests.post',  unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            self.test_data_directory,
            'Netmhcstab.other_error.html'
           ))), self.assertRaises(Exception) as context:
            output_file = tempfile.NamedTemporaryFile()
            lib.netmhc_stab.main([
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name
            ])
        self.assertTrue('Unexpected return value from NetMHCstabpan server. Unable to parse response.' in str(context.exception))
