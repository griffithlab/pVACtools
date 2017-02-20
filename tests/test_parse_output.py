import unittest
import os
import sys
import tempfile
from subprocess import call
import py_compile

def compare(path1, path2):
    r1 = open(path1)
    r2 = open(path2)
    result = not len(set(r1.readlines())^set(r2.readlines()))
    r1.close()
    r2.close()
    return result

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
        parse_output_input_iedb_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_21.ann.HLA-A*29:02.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_21.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_21.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_peptide_sequence_length_21.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_with_multiple_iedb_files(self):
        parse_output_input_iedb_files = [
            os.path.join(self.test_data_dir, "input.ann.HLA-A*29:02.9.tsv"),
            os.path.join(self.test_data_dir, "input.smm.HLA-A*29:02.9.tsv"),
            os.path.join(self.test_data_dir, "input.smmpmbec.HLA-A*29:02.9.tsv"),
        ]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "Test.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "Test_21.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            *parse_output_input_iedb_files,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
            '--top-score-metric=lowest',
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_Test_21.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_with_multiple_iedb_files_and_top_scores(self):
        parse_output_input_iedb_files = [
            os.path.join(self.test_data_dir, "input.ann.HLA-A*29:02.9.tsv"),
            os.path.join(self.test_data_dir, "input.smm.HLA-A*29:02.9.tsv"),
            os.path.join(self.test_data_dir, "input.smmpmbec.HLA-A*29:02.9.tsv"),
        ]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "Test.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "Test_21.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            *parse_output_input_iedb_files,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
            '--top-result-per-mutation',
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_Test_21.iedb.parsed.top.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_repetitive_deletion_at_beginning_of_sequence(self):
        parse_output_input_iedb_file = os.path.join(self.test_data_dir, "pat27_4.ann.HLA-A*02:01.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "pat27_4.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "pat27_4_18.fa.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_pat27_4_18.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_repetitive_insertion_at_beginning_of_sequence(self):
        parse_output_input_iedb_file = os.path.join(self.test_data_dir, "pat126.ann.HLA-A*01:01.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "pat126.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "pat126_17.fa.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_pat126_17.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_frameshift_variant_feature_elongation_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_elongation.ann.HLA-A*29:02.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_elongation.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_elongation.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_frameshift_variant_feature_elongation.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_frameshift_variant_feature_truncation_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation.ann.HLA-A*29:02.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_frameshift_variant_feature_truncation.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_frameshift_variant_feature_truncation2_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation2.ann.HLA-E*01:01.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation2.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation2.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_frameshift_variant_feature_truncation2.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_deletion_aa_deletion_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_deletion.ann.HLA-A*29:02.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_deletion.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_deletion.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_deletion_aa_deletion.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_deletion_aa_replacement_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_replacement.ann.HLA-A*29:02.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_replacement.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_replacement.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_deletion_aa_replacement.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_insertion_aa_insertion_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_insertion.ann.HLA-A*29:02.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_insertion.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_insertion.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_insertion_aa_insertion.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_insertion_aa_replacement_gets_parsed_correctly(self):
        parse_output_input_iedb_file  = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_replacement.ann.HLA-A*29:02.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_replacement.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_replacement.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_insertion_aa_replacement.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_class_ii(self):
        parse_output_input_iedb_file = os.path.join(self.test_data_dir, "input.nn_align.H2-IAb.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_31.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_31.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_nn_align.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_duplicate_transcripts(self):
        parse_output_input_iedb_file = os.path.join(self.test_data_dir, "input_multiple_transcripts_per_alt.ann.HLA-A*29:02.9.tsv")
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_multiple_transcripts_per_alt.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_multiple_transcripts_per_alt.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        self.assertFalse(call([
            self.python,
            self.executable,
            parse_output_input_iedb_file,
            parse_output_input_tsv_file,
            parse_output_key_file,
            parse_output_output_file.name,
        ], shell=False))

        expected_output_file  = os.path.join(self.test_data_dir, "output_multiple_transcripts_per_alt.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))
