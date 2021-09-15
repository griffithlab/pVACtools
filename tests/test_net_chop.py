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

from .test_utils import *
from pvactools.lib.net_chop import NetChop

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
            'net_chop.rejected.html'
        ), mode='rb')
        self.rejected = False
    else:
        reader = open(os.path.join(
            path,
            'net_chop.cterm.html'
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
            'net_chop.py'
        )
        cls.test_data_directory = os.path.join(
            pvactools_directory(),
            'tests',
            'test_data',
            'net_chop'
        )
        cls.test_fasta = os.path.join(cls.test_data_directory, 'Test_net_chop.fasta')
        cls.rejected = True

    def test_net_chop_compiles(self):
        compiled_script_path = py_compile.compile(self.script_path)
        self.assertTrue(compiled_script_path)

    def test_net_chop_cterm_runs(self):
        for method in ['cterm', '20s']:
            with patch('pvactools.lib.net_chop.requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
                data,
                files,
                self.test_data_directory,
                'net_chop.{}.html'.format(method)
                ))):
                output_file = tempfile.NamedTemporaryFile()
                NetChop(
                    os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                    self.test_fasta,
                    output_file.name,
                    method
                ).execute()
                self.assertTrue(cmp(
                    os.path.join(self.test_data_directory, 'output_{}.tsv'.format(method)),
                    output_file.name
                ))

    def test_net_chop_without_cleavage_scores(self):
        with patch('pvactools.lib.net_chop.requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            self.test_data_directory,
            'net_chop.no_cleavage.html'
            ))):
            output_file = tempfile.NamedTemporaryFile()
            NetChop(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                self.test_fasta,
                output_file.name,
                '20s'
            ).execute()
            self.assertTrue(cmp(
                os.path.join(self.test_data_directory, 'output_no_cleavage_score.tsv'),
                output_file.name
            ))

    def test_net_chop_fail(self):
        with patch('pvactools.lib.net_chop.requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            self.test_data_directory,
            'net_chop.fail.html'
           ))), self.assertRaises(Exception) as context:
            output_file = tempfile.NamedTemporaryFile()
            NetChop(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                self.test_fasta,
                output_file.name,
            ).execute()

        self.assertTrue('NetChop encountered an error during processing.' in str(context.exception))

    #This is to ensure that we catch error cases that are not explicitly handled
    def test_net_chop_other_error(self):
        with patch('pvactools.lib.net_chop.requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            self.test_data_directory,
            'net_chop.other_error.html'
           ))), self.assertRaises(Exception) as context:
            output_file = tempfile.NamedTemporaryFile()
            NetChop(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                self.test_fasta,
                output_file.name,
            ).execute()
        self.assertTrue('Unexpected return value from NetChop server. Unable to parse response.' in str(context.exception))

    def test_net_chop_rejected(self):
        logging.disable(logging.NOTSET)
        with patch('pvactools.lib.net_chop.requests.post',  unittest.mock.Mock(side_effect = lambda url, data, files=None: make_rejected_response(
            data,
           files,
           self.test_data_directory,
           self
           ))), LogCapture() as l:
            output_file = tempfile.NamedTemporaryFile()
            NetChop(
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                self.test_fasta,
                output_file.name,
            ).execute()
            self.assertTrue(cmp(
                os.path.join(self.test_data_directory, 'output_cterm.tsv'),
                output_file.name
            ))
            l.check_present(('root', 'WARNING', S("Too many jobs submitted to NetChop server. Waiting to retry.")))
