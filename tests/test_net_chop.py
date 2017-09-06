import unittest
import unittest.mock
import os
import re
import sys
import tempfile
import py_compile
from filecmp import cmp
import pvacseq.lib

def make_response(data, files, path):
    reader = open(os.path.join(
        path,
        'net_chop_%s.html'%data['method']
    ), mode='rb')
    response_obj = lambda :None
    response_obj.status_code = 200
    response_obj.content = reader.read()
    reader.close()
    return response_obj


def make_response_without_cleavage_positions(data, files, path):
    reader = open(os.path.join(
        path,
        'net_chop_2.html',
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
            'pvacseq',
            'lib',
            'net_chop.py'
        )
        cls.test_data_directory = os.path.join(
            pvac_dir,
            'tests',
            'test_data',
            'net_chop'
        )

    def test_net_chop_compiles(self):
        compiled_script_path = py_compile.compile(self.script_path)
        self.assertTrue(compiled_script_path)

    def test_net_chop_runs(self):
        self.methods=['cterm', '20s']
        self.request_mock = unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            self.test_data_directory
        ))
        pvacseq.lib.net_chop.requests.post = self.request_mock
        for method in self.methods:
            output_file = tempfile.NamedTemporaryFile()
            pvacseq.lib.net_chop.main([
                os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
                output_file.name,
                '--method',
                method
            ])
            self.assertTrue(cmp(
                os.path.join(self.test_data_directory, 'output_%s.tsv'%method),
                output_file.name
            ))

    def test_net_chop_without_cleavage_scores(self):
        self.request_mock = unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response_without_cleavage_positions(
            data,
            files,
            self.test_data_directory
        ))
        pvacseq.lib.net_chop.requests.post = self.request_mock
        output_file = tempfile.NamedTemporaryFile()
        pvacseq.lib.net_chop.main([
            os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
            output_file.name,
            '--method',
            '20s',
        ])
        self.assertTrue(cmp(
            os.path.join(self.test_data_directory, 'output_no_cleavage_score.tsv'),
            output_file.name
        ))
