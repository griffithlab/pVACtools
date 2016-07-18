import unittest
import unittest.mock
import os
from subprocess import run, call, PIPE
import re
import sys
import tempfile
import py_compile
from filecmp import cmp
from shutil import copyfile
pvac_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(pvac_dir)
import pvacseq.lib

import pdb

def make_response(method, allele, length, path):
    reader = open(os.path.join(
        path,
        'response_%s_%s_%s.tsv' % (allele, length, method)
    ), mode='r')
    response_obj = lambda :None
    response_obj.status_code = 200
    response_obj.text = reader.read()
    reader.close()
    return response_obj

def generate_call(method, allele, length, path, input_path):
    reader = open(os.path.join(
        input_path,
        "Test_21.fa"
    ), mode='r')
    text = reader.read()
    reader.close()
    return unittest.mock.call('http://tools-api.iedb.org/tools_api/mhci/', data={
        'sequence_text': ""+text,
        'method':        method,
        'allele':        allele,
        'length':        length,
    })

class PVACTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pVac_directory =  pvac_dir
        cls.test_data_directory = os.path.join(
            cls.pVac_directory,
            'tests',
            'test_data',
            'pvacseq'
        )
        cls.methods = {
            'ann': {
                'HLA-E*01:01': [9, 10]
            },
            'pickpocket': {
                'HLA-G*01:09': [9, 10],
                'HLA-E*01:01': [9, 10],
            },
        }
        cls.request_mock = unittest.mock.Mock(side_effect = lambda sequence, data: make_response(
            data['method'],
            data['allele'],
            int(data['length']),
            cls.test_data_directory
        ))
        pvacseq.lib.call_iedb.requests.post = cls.request_mock

    def test_pvacseq_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            'pvacseq',
            "pvacseq.py"
        ))
        self.assertTrue(compiled_pvac_path)

    def test_pvacseq_commands(self):
        pvac_script_path = os.path.join(
            self.pVac_directory,
            'pvacseq',
            "pvacseq.py"
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "convert_vcf",
            "generate_fasta",
            "binding_filter",
            "coverage_filter",
            "generate_fasta_key",
            "parse_output",
            "run",
            "install_vep_plugin",
            "download_example_data",
            "valid_alleles",
            "combine_parsed_outputs",
            ]:
            result = run([
                sys.executable,
                pvac_script_path,
                command,
                '-h'
            ], shell=False, stdout=PIPE)
            self.assertFalse(result.returncode)
            self.assertRegex(result.stdout.decode(), usage_search)

    def test_main_compiles(self):
        compiled_main_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            'pvacseq',
            "lib",
            "main.py"
        ))
        self.assertTrue(compiled_main_path)

    def test_pvacseq_pipeline(self):
        pvac_script_path = os.path.join(
            self.pVac_directory,
            'pvacseq',
            "pvacseq.py"
            )
        output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
        pvacseq.lib.main.main([
            os.path.join(self.test_data_directory, "input.vcf"),
            'Test',
            'HLA-G*01:09,HLA-E*01:01',
            '9,10',
            'NetMHC',
            'PickPocket',
            output_dir.name
        ])
        self.assertTrue(cmp(
            os.path.join(output_dir.name, "Test.tsv"),
            os.path.join(self.test_data_directory, "Test.tsv")
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, "Test_21.fa"),
            os.path.join(self.test_data_directory, "Test_21.fa"),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, "Test_21.key"),
            os.path.join(self.test_data_directory, "Test_21.key"),
            False
        ))
        self.assertEqual(len(self.request_mock.mock_calls), 6)
        methods = self.methods
        for method in methods.keys():
            for allele in methods[method].keys():
                for length in methods[method][allele]:
                    self.request_mock.assert_has_calls([
                        generate_call(method, allele, length, self.test_data_directory, output_dir.name)
                    ])
                    output_file   = os.path.join(output_dir.name, 'Test.%s.%s.%s.tsv' % (allele, length, method))
                    expected_file = os.path.join(self.test_data_directory, 'Test.%s.%s.%s.tsv' % (allele, length, method))
                    self.assertTrue(cmp(output_file, expected_file, False))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.HLA-E*01:01.9.parsed.tsv'),
            os.path.join(self.test_data_directory, 'Test.HLA-E*01:01.9.parsed.tsv'),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.HLA-E*01:01.10.parsed.tsv'),
            os.path.join(self.test_data_directory, 'Test.HLA-E*01:01.10.parsed.tsv'),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.HLA-G*01:09.9.parsed.tsv'),
            os.path.join(self.test_data_directory, 'Test.HLA-G*01:09.9.parsed.tsv'),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.HLA-G*01:09.10.parsed.tsv'),
            os.path.join(self.test_data_directory, 'Test.HLA-G*01:09.10.parsed.tsv'),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.combined.parsed.tsv'),
            os.path.join(self.test_data_directory, 'Test.combined.parsed.tsv'),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, "Test_filtered.tsv"),
            os.path.join(self.test_data_directory, "Test_filtered.tsv"),
            False
        ))
        output_dir.cleanup()

    def test_split_file(self):
        import random
        random.seed()
        source_file = tempfile.NamedTemporaryFile()
        writer = open(source_file.name, mode='w')
        writer.writelines(
            "".join(chr(random.randint(32,255)) for i in range(25))+"\n"
            for line in range(500)
        )
        writer.close()
        for trial in range(5):
            reader = open(source_file.name, mode='r')
            counter = 0
            total = 0
            split = random.randint(10,500)
            for chunk in pvacseq.lib.main.split_file(reader, split):
                lines = len([line for line in chunk])
                self.assertLessEqual(lines, split)
                total+=lines
            self.assertEqual(total, 500)
            reader.close()
