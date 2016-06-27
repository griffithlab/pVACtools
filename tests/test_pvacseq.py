import unittest
import os
from subprocess import run, PIPE
import re
import sys
import py_compile

class PVACTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pVac_directory = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

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
