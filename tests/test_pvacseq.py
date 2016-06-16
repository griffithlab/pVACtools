import unittest
import os
from subprocess import run, call, PIPE
import re
import sys
import tempfile
import py_compile
from filecmp import cmp
from shutil import copyfile

class PVACTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pVac_directory = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
        cls.test_data_directory = os.path.join(
            cls.pVac_directory,
            'tests',
            'test_data',
            'pvacseq'
        )

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
            "generate_fasta_key",
            "parse_output",
            "run",
            "install_vep_plugin",
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
        pvac_pipeline_cmd = "%s %s run %s Test fake_path HLA-A29:02 9 %s" % (
            sys.executable,
            pvac_script_path,
            os.path.join(self.test_data_directory, "input.vcf"),
            output_dir.name
        )
        copyfile(
            os.path.join(self.test_data_directory, 'Test.HLA-A29:02.9.netmhc.xls'),
            os.path.join(output_dir.name, 'Test.HLA-A29:02.9.netmhc.xls')
        )
        result = call([pvac_pipeline_cmd], shell=True, stdout=PIPE)
        self.assertFalse(result)
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
        self.assertTrue(cmp(
            os.path.join(output_dir.name, 'Test.HLA-A29:02.9.netmhc.parsed'),
            os.path.join(self.test_data_directory, 'Test.HLA-A29:02.9.netmhc.parsed'),
            False
        ))
        self.assertTrue(cmp(
            os.path.join(output_dir.name, "Test_filtered.xls"),
            os.path.join(self.test_data_directory, "Test_filtered.xls"),
            False
        ))
        output_dir.cleanup()
