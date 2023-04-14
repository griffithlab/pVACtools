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

from pvactools.lib.netmhc_stab import NetMHCStab
from tests.utils import *

def default_alleles(alleles):
    return ['HLA-G*01:09', 'HLA-E*01:01']

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
            "Netmhcstab.{}.html".format(data['allele'])
        ), mode='rb')
    response_obj = lambda :None
    response_obj.status_code = 200
    response_obj.content = reader.read()
    reader.close()
    return response_obj

def make_invalid_allele_response(data, files, path):
    if data['allele'] == 'HLA-B39:90':
        reader = open(os.path.join(
            path,
            "Netmhcstab.invalid_allele.html"
        ), mode='rb')
    else:
        reader = open(os.path.join(
            path,
            "Netmhcstab.{}.html".format(data['allele'])
        ), mode='rb')
    response_obj = lambda :None
    response_obj.status_code = 200
    response_obj.content = reader.read()
    reader.close()
    return response_obj

class NetChopTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.script_path = os.path.join(
            pvactools_directory(),
            'pvactools',
            'lib',
            'netmhc_stab.py'
        )
        cls.test_data_directory = os.path.join(
            pvactools_directory(),
            'tests',
            'test_data',
            'netmhc_stab'
        )
        cls.rejected = True

    def test_netmhc_stab_compiles(self):
        compiled_script_path = py_compile.compile(self.script_path)
        self.assertTrue(compiled_script_path)

    def test_netmhc_stab_runs(self):
        with patch('pvactools.lib.netmhc_stab.requests.post', unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            "Netmhcstab.{}.html".format(data['allele'])
        ))), unittest.mock.patch('pvactools.lib.netmhc_stab.NetMHCStab.valid_alleles', side_effect=default_alleles):
            output_file = tempfile.NamedTemporaryFile()
            NetMHCStab(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name,
                file_type='pVACseq'
            ).execute()
            self.assertTrue(cmp(
                os.path.join(self.test_data_directory, 'Test_filtered.stab.tsv'),
                output_file.name
            ))

    def test_netmhc_stab_fail(self):
        with patch('pvactools.lib.netmhc_stab.requests.post', unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            "Netmhcstab.fail.html"
        ))), self.assertRaises(Exception) as context, unittest.mock.patch('pvactools.lib.netmhc_stab.NetMHCStab.valid_alleles', side_effect=default_alleles):
            output_file = tempfile.NamedTemporaryFile()
            NetMHCStab(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name,
                file_type='pVACseq'
            ).execute()
        self.assertTrue('NetMHCstabpan encountered an error during processing.' in str(context.exception))

    def test_netmhc_stab_rejected(self):
        logging.disable(logging.NOTSET)
        with patch('pvactools.lib.netmhc_stab.requests.post',  unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: make_rejected_response(
            data,
            files,
            self.test_data_directory,
            self
          ))), LogCapture() as l, unittest.mock.patch('pvactools.lib.netmhc_stab.NetMHCStab.valid_alleles', side_effect=default_alleles):
            output_file = tempfile.NamedTemporaryFile()
            NetMHCStab(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name,
                file_type='pVACseq'
            ).execute()
            self.assertTrue(cmp(
                os.path.join(self.test_data_directory, 'Test_filtered.stab.tsv'),
                output_file.name
            ))
            l.check_present(('root', 'WARNING', S("Too many jobs submitted to NetMHCstabpan server. Waiting to retry.")))

    def test_netmhc_stab_cannot_open_file_error_retry(self):
        with patch('pvactools.lib.netmhc_stab.requests.post', unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            'Netmhcstab.cannot_open_file_error.html'
           ))), self.assertRaises(Exception) as context, unittest.mock.patch('pvactools.lib.netmhc_stab.NetMHCStab.valid_alleles', side_effect=default_alleles):
            output_file = tempfile.NamedTemporaryFile()
            NetMHCStab(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name,
                file_type='pVACseq'
            ).execute()
        self.assertTrue('NetMHCstabpan server was unable to read the submitted fasta file:' in str(context.exception))

    #This is to ensure that we catch error cases that are not explicitly handled
    def test_netmhc_stab_other_error(self):
        with patch('pvactools.lib.netmhc_stab.requests.post', unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            'Netmhcstab.other_error.html'
           ))), self.assertRaises(Exception) as context, unittest.mock.patch('pvactools.lib.netmhc_stab.NetMHCStab.valid_alleles', side_effect=default_alleles):
            output_file = tempfile.NamedTemporaryFile()
            NetMHCStab(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name,
                file_type='pVACseq'
            ).execute()
        self.assertTrue('Unexpected return value from NetMHCstabpan server. Unable to parse response.' in str(context.exception))

    def test_invalid_alleles(self):
        with patch('pvactools.lib.netmhc_stab.requests.post',  unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: make_invalid_allele_response(
            data,
            files,
            self.test_data_directory,
          ))):
            valid_alleles = NetMHCStab(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                "",
                file_type='pVACseq'
            ).valid_alleles(['HLA-G*01:09', 'HLA-E*01:01', 'HLA-B*39:90'])
            self.assertEqual(valid_alleles, ['HLA-G*01:09', 'HLA-E*01:01'])
