import unittest
import os
import sys
import tempfile
import py_compile

from pvactools.lib.output_parser import DefaultOutputParser, FusionOutputParser, UnmatchedSequencesOutputParser
from .test_utils import *

class OutputParserTests(unittest.TestCase):
    @classmethod
    def setUp(cls):
        executable_dir    = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable    = os.path.join(executable_dir, 'output_parser.py')
        cls.test_data_dir = os.path.join(pvactools_directory(), 'tests', 'test_data', 'output_parser')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_parse_output_runs_and_produces_expected_output(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_peptide_sequence_length_21.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_21.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_21.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_peptide_sequence_length_21',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
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

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_files,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_Test_21.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_repetitive_deletion_at_beginning_of_sequence(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "pat27_4.ann.HLA-A*02:01.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "pat27_4.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "pat27_4_18.fa.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'pat27_4',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_pat27_4_18.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_repetitive_insertion_at_beginning_of_sequence(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "pat126.ann.HLA-A*01:01.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "pat126.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "pat126_17.fa.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'pat126',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_pat126_17.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_frameshift_variant_feature_elongation_gets_parsed_correctly(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_frameshift_variant_feature_elongation.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_elongation.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_elongation.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_frameshift_variant_feature_elongation',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_frameshift_variant_feature_elongation.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_frameshift_variant_feature_truncation_gets_parsed_correctly(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_frameshift_variant_feature_truncation',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_frameshift_variant_feature_truncation.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_frameshift_variant_feature_truncation2_gets_parsed_correctly(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation2.ann.HLA-E*01:01.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation2.tsv")
        parse_output_key_file  = os.path.join(self.test_data_dir, "input_frameshift_variant_feature_truncation2.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_frameshift_variant_feature_truncation2',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())

        expected_output_file  = os.path.join(self.test_data_dir, "output_frameshift_variant_feature_truncation2.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_frameshift_variant_position_1_gets_parsed_correctly(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_frameshift_variant_position_1.MHCnuggetsI.HLA-A*02:01.8.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_frameshift_variant_position_1.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_frameshift_variant_position_1.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_frameshift_variant_position_1',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_frameshift_variant_position_1.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_deletion_aa_deletion_gets_parsed_correctly(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_inframe_deletion_aa_deletion.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_deletion.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_deletion.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_inframe_deletion_aa_deletion',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_deletion_aa_deletion.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_deletion_aa_replacement_gets_parsed_correctly(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_inframe_deletion_aa_replacement.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_replacement.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_inframe_deletion_aa_replacement.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_inframe_deletion_aa_replacement',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_deletion_aa_replacement.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_insertion_aa_insertion_gets_parsed_correctly(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_inframe_insertion_aa_insertion.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_insertion.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_insertion.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_inframe_insertion_aa_insertion',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_insertion_aa_insertion.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_input_inframe_insertion_aa_replacement_gets_parsed_correctly(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_inframe_insertion_aa_replacement.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_replacement.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_inframe_insertion_aa_replacement.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_inframe_insertion_aa_replacement',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_inframe_insertion_aa_replacement.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_class_ii(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input.nn_align.H2-IAb.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_31.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_peptide_sequence_length_31.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_nn_align.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_duplicate_transcripts(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_multiple_transcripts_per_alt.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_multiple_transcripts_per_alt.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_multiple_transcripts_per_alt.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_multiple_transcripts_per_alt',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_multiple_transcripts_per_alt.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_mnps(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_mnp.ann.HLA-A*01:01.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_mnp.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_mnp.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_mnp',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_mnp.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_mnp_at_beginning_of_sequence(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_mnp2.ann.HLA-A*01:01.10.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_mnp2.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_mnp2.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_mnp2',
        }
        parser = DefaultOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_mnp2.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_with_iedb_dna_warning(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_iedb_dna_warning.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "Test.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "Test_21.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_iedb_dna_warning',
        }
        parser = DefaultOutputParser(**parse_output_params)
        self.assertFalse(parser.execute())

    def test_parse_output_runs_and_produces_expected_output_for_fusions(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_fusions.ann.HLA-A*29:02.9.tsv")]
        parse_output_input_tsv_file = os.path.join(self.test_data_dir, "input_fusions.tsv")
        parse_output_key_file = os.path.join(self.test_data_dir, "input_fusions.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : parse_output_input_tsv_file,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_fusions',
        }
        parser = FusionOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_fusions.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_pvacvector(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_pvacvector.ann.H-2-Kb.8.tsv")]
        parse_output_key_file = os.path.join(self.test_data_dir, "input_pvacvector.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : None,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_pvacvector',
        }
        parser = UnmatchedSequencesOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_pvacvector.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_for_none_percentile(self):
        parse_output_input_iedb_file = [os.path.join(self.test_data_dir, "input_percentile_none.netmhcpan.HLA-C*03:03.9.tsv_1-2")]
        parse_output_key_file = os.path.join(self.test_data_dir, "input_pvacvector.key")
        parse_output_output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'input_iedb_files'       : parse_output_input_iedb_file,
            'input_tsv_file'         : None,
            'key_file'               : parse_output_key_file,
            'output_file'            : parse_output_output_file.name,
            'sample_name'            : 'input_percentile_none',
        }
        parser = UnmatchedSequencesOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output_percentile_none.iedb.parsed.tsv")
        self.assertTrue(compare(parse_output_output_file.name, expected_output_file))
