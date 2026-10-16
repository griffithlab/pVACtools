import unittest
import os
import sys
import tempfile
from subprocess import call
from subprocess import run as subprocess_run
from subprocess import PIPE
from filecmp import cmp
import pandas as pd
import py_compile

from tests.utils import *
from pvactools.lib.color_peptides51mer import annotate_every_nucleotide, set_underline, get_mutant_positions_from_fasta


class PvacspliceCreatePeptideOrderingFormTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(
            pvactools_directory(), "pvactools", "tools", "pvacsplice"
        )
        cls.executable = os.path.join(
            cls.executable_dir, "create_peptide_ordering_form.py"
        )
        cls.test_data_dir = os.path.join(
            pvactools_directory(),
            "tests",
            "test_data",
            "pvacsplice_create_peptide_ordering_form",
        )

        cls.transcripts_fasta = os.path.join(cls.test_data_dir, "HCC1395_TUMOR_DNA.transcripts.fa")
        cls.mhc_class_i_tsv = os.path.join(
            cls.test_data_dir, "MHC_Class_I", "HCC1395_TUMOR_DNA.MHC_I.all_epitopes.aggregated.tsv"
        )
        cls.mhc_class_ii_tsv = os.path.join(
            cls.test_data_dir, "MHC_Class_II", "HCC1395_TUMOR_DNA.MHC_II.all_epitopes.aggregated.tsv"
        )
        cls.sample_name = "HCC1395_TUMOR_DNA"

        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_file_prefix = "output"
        cls.output_path = cls.tmpdir.name

        result = call(
            [
                cls.python,
                cls.executable,
                cls.transcripts_fasta,
                '25',
                cls.mhc_class_i_tsv,
                cls.mhc_class_ii_tsv,
                cls.output_file_prefix,
                cls.sample_name,
                "-o",
                cls.output_path,
                "--aggregate-report-evaluation",
                "Pending",
            ],
            shell=False,
        )

        if result != 0:
            cls.tmpdir.cleanup()
            raise RuntimeError("Failed to run create_peptide_ordering_form.py")

    @classmethod
    def tearDownClass(cls):
        cls.tmpdir.cleanup()

    def test_command(self):
        pvac_script_path = os.path.join(
            self.executable_dir,
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        result = subprocess_run([
            sys.executable,
            pvac_script_path,
            'create_peptide_ordering_form',
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacsplice create_peptide_ordering_form -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

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

    def test_frameshift_splice_site_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "RDFCFGPWKLTASKTHIMKSADVVKQRFKNPAWVWLWN",
            classI_peptide = "RFKNPAWVWL",
            classII_peptide = "HIMKSADVVKQRFKN",
            classI_ic50 = "5620.54",
            classI_percentile = "1.8",
            classII_ic50 = "742.4",
            classII_percentile = "30.0",
            classI_transcript = "ENST00000367833.7",
            classII_transcript = "ENST00000367833.7",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_frameshift_splice_site.fa")
        full_id = "21.TIPRL.ENST00000367833.JUNC00000028.chr1:168179177-168179178.D.frameshift_splice_site"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {32, 33, 34, 35, 36, 37, 25, 26, 27, 28, 29, 30, 31}

        # Assert all mutant AAs are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )

        # Check color applied to class I peptide
        color_indices = [26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
        for idx in color_indices:
            assert peptide_sequence[idx].color, f"Amino acid at index {idx} should be colored"


    def test_inframe_splice_site_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "LVLHLNQLEGNKEKFEKQLKKKSEEKELKIKNHSLQETSEQNVILQHTLQ",
            classI_peptide = "KSEEKELKI",
            classII_peptide = "KFEKQLKKKSEEKEL",
            classI_ic50 = "9759.45",
            classI_percentile = "20.0",
            classII_ic50 = "56.1",
            classII_percentile = "5.4",
            classI_transcript = "ENST00000690025.1",
            classII_transcript = "ENST00000690025.1",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_inframe_splice_site.fa")
        full_id = "15.CCDC18.ENST00000690025.JUNC00000011.chr1:93226448-93226449.D.inframe_splice_site"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {24, 25}

        # Assert all mutant AAs are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )

        # Check bold applied to class II peptide
        bold_indices = [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]
        for idx in bold_indices:
            assert peptide_sequence[idx].bold, f"Amino acid at index {idx} should be bolded"
