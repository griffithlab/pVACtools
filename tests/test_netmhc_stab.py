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
            'pvacseq',
            'lib',
            'netmhc_stab.py'
        )
        cls.test_data_directory = os.path.join(
            pvac_dir,
            'tests',
            'test_data',
            'netmhc_stab'
        )
        cls.request_mock = unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            cls.test_data_directory
        ))
        pvacseq.lib.netmhc_stab.requests.post = cls.request_mock

    def test_netmhc_stab_compiles(self):
        compiled_script_path = py_compile.compile(self.script_path)
        self.assertTrue(compiled_script_path)

    def test_netmhc_stab_runs(self):
        output_file = tempfile.NamedTemporaryFile()
        pvacseq.lib.netmhc_stab.main([
            os.path.join(self.test_data_directory, 'Test_filtered.tsv'),
            output_file.name
        ])
        self.assertTrue(cmp(
            os.path.join(self.test_data_directory, 'Test_filtered.stab.tsv'),
            output_file.name
        ))
