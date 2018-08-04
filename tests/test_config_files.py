import unittest
import sys
import os
import py_compile
from subprocess import check_output

class ConfigFilesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'tools', 'pvacseq')
        cls.executable     = os.path.join(cls.executable_dir, 'config_files.py')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_pvacseq_config_files_additional_input_file_list_command(self):
        output = check_output(
            [self.python, self.executable, 'additional_input_file_list'],
            shell = False
        )
        expected_output = (
            "gene_expn_file: <genes.fpkm_tracking file from Cufflinks>\n"
            + "transcript_expn_file: <isoforms.fpkm_tracking file from Cufflinks>\n"
            + "normal_snvs_coverage_file: <bam-readcount output file for normal BAM and snvs>\n"
            + "normal_indels_coverage_file: <bam-readcount output file for normal BAM and indels>\n"
            + "tdna_snvs_coverage_file: <bam-readcount output file for tumor DNA BAM and snvs>\n"
            + "tdna_indels_coverage_file: <bam-readcount output file for tumor DNA BAM and indels>\n"
            + "trna_snvs_coverage_file: <bam-readcount output file for tumor RNA BAM and snvs>\n"
            + "trna_indels_coverage_file: <bam-readcount output file for tumor RNA BAM and indels>\n"
            + "phased_proximal_variants_vcf: <A VCF with phased proximal variant information>\n"
        )
        self.assertEqual(output.decode(), expected_output)
