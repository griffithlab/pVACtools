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

def make_response(method, path):
    reader = open(os.path.join(
        path,
        'response_%s.tsv' %method
    ), mode='r')
    response_obj = lambda :None
    response_obj.status_code = 200
    response_obj.text = reader.read()
    reader.close()
    return response_obj

def generate_call(method, path, input_path):
    reader = open(os.path.join(
        input_path,
        "Test_21.fa"
    ), mode='r')
    text = reader.read()
    reader.close()
    return unittest.mock.call('http://tools-api.iedb.org/tools_api/mhci/', data={
        'sequence_text': ""+text,
        'method':        method,
        'allele':        "HLA-A*29:02",
        'length':        9,
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
        cls.methods = ['ann', 'smm', 'smmpmbec']
        cls.request_mock = unittest.mock.Mock(side_effect = (
            make_response(method, cls.test_data_directory) for method in cls.methods
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
            ]:
            temp_cmd = "%s %s %s -h" %(
                sys.executable,
                pvac_script_path,
                command
            )
            result = run([temp_cmd], shell=True, stdout=PIPE)
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
        pvac_pipeline_cmd = "%s %s run %s Test HLA-A29:02 9 NetMHC SMM SMMPMBEC %s" % (
            sys.executable,
            pvac_script_path,
            os.path.join(self.test_data_directory, "input.vcf"),
            output_dir.name
        )
        pvacseq.lib.main.main([
            os.path.join(self.test_data_directory, "input.vcf"),
            'Test',
            'HLA-A29:02',
            '9',
            'NetMHC',
            'SMM',
            'SMMPMBEC',
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
        self.assertEqual(len(self.request_mock.mock_calls), 3)
        self.request_mock.assert_has_calls([
            generate_call(method, self.test_data_directory, output_dir.name)
            for method in self.methods
        ])
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.HLA-A29:02.9.ann.tsv'),
            os.path.join(self.test_data_directory, 'Test.HLA-A29:02.9.ann.tsv'),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.HLA-A29:02.9.smm.tsv'),
            os.path.join(self.test_data_directory, 'Test.HLA-A29:02.9.smm.tsv'),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.HLA-A29:02.9.smmpmbec.tsv'),
            os.path.join(self.test_data_directory, 'Test.HLA-A29:02.9.smmpmbec.tsv'),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.HLA-A29:02.9.parsed.tsv'),
            os.path.join(self.test_data_directory, 'Test.HLA-A29:02.9.parsed.tsv'),
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
        test_dir = tempfile.TemporaryDirectory()
        prefix = os.path.join(test_dir.name, 'split_file_')
        writer = open(prefix+"source", mode='w')
        writer.writelines(
            "".join(chr(random.randint(32,255)) for i in range(25))+"\n"
            for line in range(500)
        )
        writer.close()
        for trial in range(5):
            reader = open(prefix+"source", mode='r')
            counter = 0
            split = random.randint(10,500)
            for chunk in pvacseq.lib.main.split_file(reader, split):
                writer = open(prefix+"output_%d_%d"%(trial, counter), mode='w')
                counter+=1
                writer.writelines(chunk)
                writer.close()
            total = 0
            reader.close()
            for i in range(counter):
                reader = open(prefix+"output_%d_%d"%(trial, i), mode='r')
                lines = len(reader.read().strip().split("\n"))
                reader.close()
                self.assertLessEqual(lines, split)
                total+=lines
            self.assertEqual(total, 500)
        test_dir.cleanup()
