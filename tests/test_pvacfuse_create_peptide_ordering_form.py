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


class PvacfuseCreatePeptideOrderingFormTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(
            pvactools_directory(), "pvactools", "tools", "pvacfuse"
        )
        cls.executable = os.path.join(
            cls.executable_dir, "create_peptide_ordering_form.py"
        )
        cls.test_data_dir = os.path.join(
            pvactools_directory(),
            "tests",
            "test_data",
            "pvacfuse_create_peptide_ordering_form",
        )

        cls.input = os.path.join(cls.test_data_dir, "agfusion_HCC1395")
        cls.peptide_fasta = os.path.join(cls.test_data_dir, "Homo_sapiens.GRCh38.pep.short.fa.gz")
        transcript_fasta = os.path.join(cls.test_data_dir, 'Homo_sapiens.GRCh38.95.cds.all.fa.gz')
        cls.unzipped_transcript_fasta = gunzip_file(transcript_fasta)
        cls.mhc_class_i_tsv = os.path.join(
            cls.test_data_dir, "MHC_Class_I", "sample.name.MHC_I.all_epitopes.aggregated.tsv"
        )
        cls.mhc_class_ii_tsv = os.path.join(
            cls.test_data_dir, "MHC_Class_II", "sample.name.MHC_II.all_epitopes.aggregated.tsv"
        )
        cls.sample_name = "HCC1395_TUMOR_DNA"

        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_file_prefix = "output"
        cls.output_path = cls.tmpdir.name

        result = call(
            [
                cls.python,
                cls.executable,
                cls.input,
                cls.unzipped_transcript_fasta,
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
        self.assertFalse(result.returncode, "Failed `pvacfuse create_peptide_ordering_form -h`")
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

    def test_frameshift_fusion_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "QEERPIRQILYLGDLLETCHFQAFWHQFQRMSF",
            classI_peptide = "HFQAFWHQF",
            classII_peptide = "CHFQAFWHQFQRMSF",
            classI_ic50 = "257.23",
            classI_percentile = "0.73",
            classII_ic50 = "145.1",
            classII_percentile = "12.0",
            classI_transcript = "ENST00000588934-ENST00000275016",
            classII_transcript = "ENST00000588934-ENST00000275016",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_frameshift_fusion.fa")
        full_id = "7.EIF3K_CYP39A1.ENST00000588934_ENST00000275016.frameshift_fusion.118"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {32, 25, 26, 27, 28, 29, 30, 31}

        # Assert all mutant AAs are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )

        # Check color applied to class I peptide
        color_indices = [19, 20, 21, 22, 23, 24, 25, 26, 27]
        for idx in color_indices:
            assert peptide_sequence[idx].color, f"Amino acid at index {idx} should be colored"

        # Check bold applied to class II peptide
        bold_indices = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
        for idx in bold_indices:
            assert peptide_sequence[idx].bold, f"Amino acid at index {idx} should be bolded"


    def test_inframe_fusion_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "ADTKYNDTDRCWGPLRRVDAYRIYLTILVGDSGVGKTSLLVQFDQGKFIP",
            classI_peptide = "YRIYLTILV",
            classII_peptide = "PLRRVDAYRIYLTIL",
            classI_ic50 = "6204.38",
            classI_percentile = "5.0",
            classII_ic50 = "292.6",
            classII_percentile = "19.0",
            classI_transcript = "ENST00000417024-ENST00000340415",
            classII_transcript = "ENST00000417024-ENST00000340415",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_inframe_fusion.fa")
        full_id = "23.TMEM104_RAB37.ENST00000417024_ENST00000340415.inframe_fusion.225"

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
