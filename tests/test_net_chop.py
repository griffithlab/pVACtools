import unittest
import unittest.mock
import os
import re
import sys
import tempfile
import py_compile
from filecmp import cmp
pvac_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(pvac_dir)
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

class NetChopTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
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
        cls.methods=['cterm', '20s']
        cls.request_mock = unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            cls.test_data_directory
        ))
        pvacseq.lib.net_chop.requests.post = cls.request_mock

    def test_net_chop_compiles(self):
        compiled_script_path = py_compile.compile(self.script_path)
        self.assertTrue(compiled_script_path)

    def test_net_chop_runs(self):
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
