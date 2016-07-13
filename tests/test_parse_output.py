import unittest
import os
import sys
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile

class ParseOutputTests(unittest.TestCase):
    @classmethod
    def setUp(cls):
        cls.python        = sys.executable
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        executable_dir    = os.path.join(base_dir, 'pvacseq', 'lib')
        cls.executable    = os.path.join(executable_dir, 'parse_output.py')
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data', 'parse_output')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_parse_output_runs_and_produces_expected_output(self):
        parse_output_input_iedb_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_21.HLA-A29:02.9.ann.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_21.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_21.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, parse_output_input_iedb_file, parse_output_input_tsv_file, parse_output_key_file, parse_output_output_file.name
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_peptide_sequence_length_21.iedb.parsed.tsv")
        self.assertTrue(cmp(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_with_multiple_iedb_files(self):
        parse_output_input_iedb_files = [
            os.path.join(self.test_data_dir, "input.HLA-A29:02.9.ann.tsv"),
            os.path.join(self.test_data_dir, "input.HLA-A29:02.9.smm.tsv"),
            os.path.join(self.test_data_dir, "input.HLA-A29:02.9.smmpmbec.tsv"),
        ]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "Test.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "Test_21.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, *parse_output_input_iedb_files, parse_output_input_tsv_file, parse_output_key_file, parse_output_output_file.name
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_Test_21.iedb.parsed.tsv")
        self.assertTrue(cmp(parse_output_output_file.name, expected_output_file))

    def test_input_frameshift_variant_feature_elongation_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_elongation.HLA-A29:02.9.ann.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_elongation.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_elongation.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, parse_output_input_iedb_file, parse_output_input_tsv_file, parse_output_key_file, parse_output_output_file.name
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_frameshift_variant_feature_elongation.iedb.parsed.tsv")
        self.assertTrue(cmp(parse_output_output_file.name, expected_output_file))

    def test_input_frameshift_variant_feature_truncation_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation.HLA-A29:02.9.ann.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, parse_output_input_iedb_file, parse_output_input_tsv_file, parse_output_key_file, parse_output_output_file.name
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_frameshift_variant_feature_truncation.iedb.parsed.tsv")
        self.assertTrue(cmp(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_deletion_aa_deletion_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_deletion.HLA-A29:02.9.ann.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_deletion.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_deletion.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, parse_output_input_iedb_file, parse_output_input_tsv_file, parse_output_key_file, parse_output_output_file.name
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_deletion_aa_deletion.iedb.parsed.tsv")
        self.assertTrue(cmp(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_deletion_aa_replacement_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_replacement.HLA-A29:02.9.ann.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_replacement.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_replacement.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, parse_output_input_iedb_file, parse_output_input_tsv_file, parse_output_key_file, parse_output_output_file.name
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_deletion_aa_replacement.iedb.parsed.tsv")
        self.assertTrue(cmp(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_insertion_aa_insertion_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_insertion.HLA-A29:02.9.ann.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_insertion.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_insertion.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, parse_output_input_iedb_file, parse_output_input_tsv_file, parse_output_key_file, parse_output_output_file.name
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_insertion_aa_insertion.iedb.parsed.tsv")
        self.assertTrue(cmp(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_insertion_aa_replacement_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_replacement.HLA-A29:02.9.ann.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_replacement.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_replacement.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python, self.executable, parse_output_input_iedb_file, parse_output_input_tsv_file, parse_output_key_file, parse_output_output_file.name
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_insertion_aa_replacement.iedb.parsed.tsv")
        self.assertTrue(cmp(parse_output_output_file.name, expected_output_file))
