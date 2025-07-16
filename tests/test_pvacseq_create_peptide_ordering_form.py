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
        self.assertFalse(result.returncode, "Failed `pvacseq create_peptide_ordering_form -h`")
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
    
    def test_annotate_classi_and_classii_pass(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence="ABCDEFGHIJKL",
            classI_peptide="DEFG",
            classII_peptide="GHI",
            classI_ic50="100",
            classI_percentile="0.4",
            classII_ic50="40",
            classII_percentile="0.2",
            classI_transcript="T1",
            classII_transcript="T1",
            cIIC50_threshold=1000,
            cIpercentile_threshold=2,
            cIIIC50_threshold=500,
            cIIpercent_threshold=2,
            probPos=""
        )

        # Check color applied to class I peptide ("DEFG")
        color_indices = [3, 4, 5, 6]
        for idx in color_indices:
            assert peptide_sequence[idx].color, f"Amino acid at index {idx} should be colored"

        # Check bold applied to class II peptide ("GHI")
        bold_indices = [6, 7, 8]
        for idx in bold_indices:
            assert peptide_sequence[idx].bold, f"Amino acid at index {idx} should be bolded"
    
    def test_annotate_classi_pass(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence="ABCDEFGHIJKL",
            classI_peptide="DEFG",
            classII_peptide="GHI",
            classI_ic50="100",
            classI_percentile="0.4",
            classII_ic50="1500",
            classII_percentile="2.2",
            classI_transcript="T1",
            classII_transcript="T1",
            cIIC50_threshold=1000,
            cIpercentile_threshold=2,
            cIIIC50_threshold=500,
            cIIpercent_threshold=2,
            probPos=""
        )

        # Check color applied to class I peptide ("DEFG")
        color_indices = [3, 4, 5, 6]
        for idx in color_indices:
            assert peptide_sequence[idx].color, f"Amino acid at index {idx} should be colored"

        # Check bold is NOT applied to class II peptide ("GHI")
        bold_indices = [6, 7, 8]
        for idx in bold_indices:
            assert not peptide_sequence[idx].bold, f"Amino acid at index {idx} should not be bolded"

    def test_annotate_classii_pass(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence="ABCDEFGHIJKL",
            classI_peptide="DEFG",
            classII_peptide="GHI",
            classI_ic50="1200",
            classI_percentile="3.0",
            classII_ic50="40",
            classII_percentile="0.2",
            classI_transcript="T1",
            classII_transcript="T1",
            cIIC50_threshold=1000,
            cIpercentile_threshold=2,
            cIIIC50_threshold=500,
            cIIpercent_threshold=2,
            probPos=""
        )

        # Check color is NOT applied to class I peptide ("DEFG")
        color_indices = [3, 4, 5, 6]
        for idx in color_indices:
            assert not peptide_sequence[idx].color, f"Amino acid at index {idx} should not be colored"

        # Check bold applied to class II peptide ("GHI")
        bold_indices = [6, 7, 8]
        for idx in bold_indices:
            assert peptide_sequence[idx].bold, f"Amino acid at index {idx} should be bolded"
    
    def test_annotate_classi_and_classii_fail(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence="ABCDEFGHIJKL",
            classI_peptide="DEFG",
            classII_peptide="GHI",
            classI_ic50="1800",
            classI_percentile="2.3",
            classII_ic50="700",
            classII_percentile="4.2",
            classI_transcript="T1",
            classII_transcript="T1",
            cIIC50_threshold=1000,
            cIpercentile_threshold=2,
            cIIIC50_threshold=500,
            cIIpercent_threshold=2,
            probPos=["C", "H"]
        )

        # Check color is NOT applied to class I peptide ("DEFG")
        color_indices = [3, 4, 5, 6]
        for idx in color_indices:
            assert not peptide_sequence[idx].color, f"Amino acid at index {idx} should not be colored"

        # Check bold is NOT applied to class II peptide ("GHI")
        bold_indices = [6, 7, 8]
        for idx in bold_indices:
            assert not peptide_sequence[idx].bold, f"Amino acid at index {idx} should not be bolded"
        
        # Check 'large' from probPos
        assert peptide_sequence[2].large  # "C"
        assert peptide_sequence[7].large  # "H"
    
    def test_missense_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "EDAVQGIANQDAAQGIAKE",
            classI_peptide = "ANQDAAQGI",
            classII_peptide = "VQGIANQDAAQGIAK",
            classI_ic50 = "800",
            classI_percentile = "1.5",
            classII_ic50 = "400",
            classII_percentile = "0.7",
            classI_transcript = "ENST00000434783.3",
            classII_transcript = "ENST00000434783.3",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_missense.fa")
        full_id = "8.FAM230A.ENST00000434783.3.missense.322E/Q"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {9}

        # Assert all mutant AAs are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )
    
    def test_missense_with_proximal_variant_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "EDASQGIANQDAAQGIAKE",
            classI_peptide = "ANQDAAQGI",
            classII_peptide = "VQGIANQDAAQGIAK",
            classI_ic50 = "800",
            classI_percentile = "1.5",
            classII_ic50 = "400",
            classII_percentile = "0.7",
            classI_transcript = "ENST00000434783.3",
            classII_transcript = "ENST00000434783.3",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_missense_proximal_variant.fa")
        full_id = "8.FAM230A.ENST00000434783.3.missense.322E/Q"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {3, 9}

        # Assert all mutant AAs are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )
    
    def test_frameshift_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "DTGGGGRSAGSTGQGSGEKAGCPWSGTGQH",
            classI_peptide = "GEKAGCPWS",
            classII_peptide = "SGEKAGCPWSGTGQH",
            classI_ic50 = "800",
            classI_percentile = "1.5",
            classII_ic50 = "400",
            classII_percentile = "0.7",
            classI_transcript = "ENST00000406386.3",
            classII_transcript = "ENST00000406386.3",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_frameshift.fa")
        full_id = "16.TRIOBP.ENST00000406386.3.FS.219GA/G"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26}

        # Assert all mutant AAs are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )
    
    def test_frameshift_with_proximal_variant_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "DTSGGGRSAGSTGQGSGEKAGCPWSGTGQH",
            classI_peptide = "GEKAGCPWS",
            classII_peptide = "SGEKAGCPWSGTGQH",
            classI_ic50 = "800",
            classI_percentile = "1.5",
            classII_ic50 = "400",
            classII_percentile = "0.7",
            classI_transcript = "ENST00000406386.3",
            classII_transcript = "ENST00000406386.3",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_frameshift_proximal_variant.fa")
        full_id = "16.TRIOBP.ENST00000406386.3.FS.219GA/G"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {2, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26}

        # Assert all mutant AAs are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )
    
    def test_inframe_deletion_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "PASAAAAAAAVIPTVSTPP",
            classI_peptide = "SAAAAAAAV",
            classII_peptide = "SAAAAAAAVIPTVST",
            classI_ic50 = "800",
            classI_percentile = "1.5",
            classII_ic50 = "400",
            classII_percentile = "0.7",
            classI_transcript = "ENST00000381793.2",
            classII_transcript = "ENST00000381793.2",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_inframe_deletion.fa")
        full_id = "2.RBM47.ENST00000381793.2.inframe_del.495-502AAAAAAAA/A"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {2, 3}

        # Assert the two AAs surrounding the deletion are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )
    
    def test_inframe_deletion_with_proximal_variant_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "PASAAAAAAAVIPTVSTPL",
            classI_peptide = "SAAAAAAAV",
            classII_peptide = "SAAAAAAAVIPTVST",
            classI_ic50 = "800",
            classI_percentile = "1.5",
            classII_ic50 = "400",
            classII_percentile = "0.7",
            classI_transcript = "ENST00000381793.2",
            classII_transcript = "ENST00000381793.2",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_inframe_deletion_proximal_variant.fa")
        full_id = "2.RBM47.ENST00000381793.2.inframe_del.495-502AAAAAAAA/A"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {2, 3, 18}

        # Assert the two AAs surrounding the deletion are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )

    def test_inframe_insertion_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "PLPPPPLLPLLPLLLLLGASGGG",
            classI_peptide = "LLPLLPLLL",
            classII_peptide = "LPLLPLLLLLGASGG",
            classI_ic50 = "800",
            classI_percentile = "1.5",
            classII_ic50 = "400",
            classII_percentile = "0.7",
            classI_transcript = "ENST00000233809.4",
            classII_transcript = "ENST00000233809.4",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_inframe_insertion.fa")
        full_id = "1.IGFBP2.ENST00000233809.4.inframe_ins.20L/LLP"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {10, 11}

        # Assert the differing inserted AAs are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )
    
    def test_inframe_insertion_with_proximal_variant_formatting(self):
        peptide_sequence = annotate_every_nucleotide(
            sequence = "PTPPPPLLPLLPLLLLLGASGGG",
            classI_peptide = "LLPLLPLLL",
            classII_peptide = "LPLLPLLLLLGASGG",
            classI_ic50 = "800",
            classI_percentile = "1.5",
            classII_ic50 = "400",
            classII_percentile = "0.7",
            classI_transcript = "ENST00000233809.4",
            classII_transcript = "ENST00000233809.4",
            cIIC50_threshold = 1000,
            cIpercentile_threshold = 2,
            cIIIC50_threshold = 500,
            cIIpercent_threshold = 2,
            probPos = ""
        )

        fasta_file = os.path.join(self.test_data_dir, "test_inframe_insertion_proximal_variant.fa")
        full_id = "1.IGFBP2.ENST00000233809.4.inframe_ins.20L/LLP"

        mutant_positions = get_mutant_positions_from_fasta(fasta_file, full_id)
        set_underline(peptide_sequence, mutant_positions)

        expected_underlined_positions = {1, 10, 11}

        # Assert the differing inserted AAs are underlined
        self.assertSetEqual(
            mutant_positions,
            expected_underlined_positions,
            f"Expected underlining at {mutant_positions}, but got {expected_underlined_positions}",
        )

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
