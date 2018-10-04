import unittest
import unittest.mock
import os
import re
import sys
import tempfile
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run
from filecmp import cmp
import yaml
import lib
import datetime
from tools.pvacfuse import *
from mock import patch
from .test_utils import *

def make_response(data, files, path):
    if not files:
        if 'length' in data:
            filename = 'response_%s_%s_%s.tsv' % (data['allele'], data['length'], data['method'])
        else:
            filename = 'response_%s_%s.tsv' % (data['allele'], data['method'])
        reader = open(os.path.join(
            path,
            filename
        ), mode='r')
        response_obj = lambda :None
        response_obj.status_code = 200
        response_obj.text = reader.read()
        reader.close()
        return response_obj
    else:
        basefile = os.path.basename(data['configfile'])
        reader = open(os.path.join(
            path,
            'net_chop.html' if basefile == 'NetChop.cf' else 'Netmhcstab.html'
        ), mode='rb')
        response_obj = lambda :None
        response_obj.status_code = 200
        response_obj.content = reader.read()
        reader.close()
        return response_obj

def generate_class_i_call(method, allele, length, input_file):
    reader = open(input_file, mode='r')
    text = reader.read()
    reader.close()
    return unittest.mock.call('http://tools-cluster-interface.iedb.org/tools_api/mhci/', data={
        'sequence_text': ""+text,
        'method':        method,
        'allele':        allele,
        'length':        length,
        'user_tool':     'pVac-seq',
    })

def generate_class_ii_call(method, allele, path, input_path):
    reader = open(os.path.join(
        input_path,
        "MHC_Class_II",
        "tmp",
        "Test_31.fa.split_1-48"
    ), mode='r')
    text = reader.read()
    reader.close()
    return unittest.mock.call('http://tools-cluster-interface.iedb.org/tools_api/mhcii/', data={
        'sequence_text': ""+text,
        'method':        method,
        'allele':        allele,
        'user_tool':     'pVac-seq',
    })

def pvac_directory():
    return os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

def test_data_directory():
    return os.path.join(
        pvac_directory(),
        'tests',
        'test_data',
        'pvacfuse'
    )

class PvacfuseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pVac_directory = pvac_directory()
        cls.test_data_directory = test_data_directory()

    def test_pvacfuse_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            'tools',
            'pvacfuse',
            'main.py'
        ))
        self.assertTrue(compiled_pvac_path)

    def test_pvacfuse_commands(self):
        pvac_script_path = os.path.join(
            self.pVac_directory,
            'tools',
            'pvacfuse',
            'main.py'
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "run",
            "binding_filter",
            "valid_alleles",
            "download_example_data",
            ]:
            result = subprocess_run([
                sys.executable,
                pvac_script_path,
                command,
                '-h'
            ], shell=False, stdout=PIPE)
            self.assertFalse(result.returncode)
            self.assertRegex(result.stdout.decode(), usage_search)

    def test_run_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            "tools",
            "pvacfuse",
            "run.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_pvacfuse_pipeline(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "fusions_annotated.bedpe"),
                'Test',
                'HLA-A*29:02',
                'NetMHC',
                output_dir.name,
                '-e', '9',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
            ])

            for file_name in (
                'Test.tsv',
                'Test.tsv_1-4',
                'Test.all_epitopes.tsv',
                'Test.filtered.tsv',
                'Test.filtered.condensed.ranked.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', file_name)
                self.assertTrue(compare(output_file, expected_file))

            for file_name in (
                'Test_21.fa.split_1-8',
                'Test_21.fa.split_1-8.key',
                'Test.ann.HLA-A*29:02.9.tsv_1-8',
                'Test.HLA-A*29:02.9.parsed.tsv_1-8',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'fusions', 'MHC_Class_I', 'tmp', file_name)
                self.assertTrue(compare(output_file, expected_file))

            mock_request.assert_has_calls([
                generate_class_i_call('ann', 'HLA-A*29:02', 9, os.path.join(output_dir.name, "MHC_Class_I", "tmp", "Test_21.fa.split_1-8"))
            ])

            output_dir.cleanup()

    def test_pvacfuse_combine_and_condense_steps(self):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
            for subdir in ['MHC_Class_I', 'MHC_Class_II']:
                path = os.path.join(output_dir.name, subdir)
                os.mkdir(path)
                test_data_dir = os.path.join(self.test_data_directory, 'combine_and_condense', subdir)
                for item in os.listdir(test_data_dir):
                    os.symlink(os.path.join(test_data_dir, item), os.path.join(path, item))

            run.main([
                os.path.join(self.test_data_directory, "fusions_annotated.bedpe"),
                'Test',
                'HLA-A*29:02,H2-IAb',
                'NetMHC', 'NNalign',
                output_dir.name,
                '-e', '9',
                '--top-score-metric=lowest',
                '--keep-tmp-files',
            ])

            for file_name in (
                'Test.all_epitopes.tsv',
                'Test.filtered.tsv',
                'Test.filtered.condensed.ranked.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'combined', file_name)
                expected_file = os.path.join(self.test_data_directory, 'combine_and_condense', 'combined', file_name)
                self.assertTrue(compare(output_file, expected_file))
