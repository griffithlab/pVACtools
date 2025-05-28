import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import pandas as pd
import py_compile

from tests.utils import *


class CreatePeptideOrderingFormTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(
            pvactools_directory(), "pvactools", "tools", "pvacseq"
        )
        cls.executable = os.path.join(
            cls.executable_dir, "create_peptide_ordering_form.py"
        )
        cls.test_data_dir = os.path.join(
            pvactools_directory(),
            "tests",
            "test_data",
            "pvacseq_create_peptide_ordering_form",
        )

        cls.input_vcf = os.path.join(cls.test_data_dir, "input.vcf")
        cls.mhc_class_i_tsv = os.path.join(
            cls.test_data_dir, "MHC_Class_I", "Test.all_epitopes.aggregated.tsv"
        )
        cls.mhc_class_ii_tsv = os.path.join(
            cls.test_data_dir, "MHC_Class_II", "Test.all_epitopes.aggregated.tsv"
        )
        cls.sample_name = "H_NJ-HCC1395-HCC1395"

        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_file_prefix = "output"
        cls.output_path = cls.tmpdir.name

        result = call(
            [
                cls.python,
                cls.executable,
                cls.input_vcf,
                cls.mhc_class_i_tsv,
                cls.mhc_class_i_tsv,
                cls.mhc_class_ii_tsv,
                cls.output_file_prefix,
                cls.sample_name,
                "-o",
                cls.output_path,
            ],
            shell=False,
        )

        if result != 0:
            cls.tmpdir.cleanup()
            raise RuntimeError("Failed to run create_peptide_ordering_form.py")

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_colored_peptides_excel_output(self):
        generated_xlsx = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}_{self.sample_name}.Colored_Peptides.xlsx",
        )
        expected_xlsx = os.path.join(
            self.test_data_dir, "expected_output.Colored_Peptides.xlsx"
        )

        self.assertTrue(os.path.exists(generated_xlsx))

        generated_df = pd.read_excel(generated_xlsx)
        expected_df = pd.read_excel(expected_xlsx)

        generated_df = (
            generated_df.sort_index(axis=1)
            .sort_values(by=list(generated_df.columns))
            .reset_index(drop=True)
        )
        expected_df = (
            expected_df.sort_index(axis=1)
            .sort_values(by=list(expected_df.columns))
            .reset_index(drop=True)
        )

        try:
            pd.testing.assert_frame_equal(generated_df, expected_df, check_dtype=False)
        except AssertionError as e:
            self.fail(f"Generated Excel content does not match expected:\n{e}")

    def test_fasta_output(self):
        generated_fasta = os.path.join(
            self.output_path, f"{self.output_file_prefix}_{self.sample_name}.fa"
        )
        expected_fasta = os.path.join(self.test_data_dir, "expected_output.fa")

        self.assertTrue(os.path.exists(generated_fasta))
        self.assertTrue(
            cmp(generated_fasta, expected_fasta), "FASTA output does not match expected"
        )

    def test_manufacturability_tsv_output(self):
        generated_tsv = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}_{self.sample_name}.manufacturability.tsv",
        )
        expected_tsv = os.path.join(
            self.test_data_dir, "expected_output.manufacturability.tsv"
        )

        self.assertTrue(os.path.exists(generated_tsv))
        self.assertTrue(
            cmp(generated_tsv, expected_tsv),
            "Manufacturability TSV does not match expected",
        )
