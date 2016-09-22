import unittest
import os
import sys
import tempfile
from filecmp import cmp
import py_compile
from pvacseq.lib.input_file_converter import *

class InputFileConverterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable_dir = os.path.join(base_dir, 'pvacseq', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'input_file_converter.py')
        cls.test_data_dir  = os.path.join(base_dir, 'tests', 'test_data', 'input_file_converter')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_input_vcf_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_cufflinks_files_generates_expected_tsv(self):
        convert_vcf_input_file              = os.path.join(self.test_data_dir, 'full_input.vcf')
        convert_vcf_output_file             = tempfile.NamedTemporaryFile()
        convert_vcf_cufflinks_genes_file    = os.path.join(self.test_data_dir, 'genes.fpkm_tracking')
        convert_vcf_cufflinks_isoforms_file = os.path.join(self.test_data_dir, 'isoforms.fpkm_tracking')

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : convert_vcf_cufflinks_genes_file,
            'transcript_expn_file'       : convert_vcf_cufflinks_isoforms_file,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_cufflinks.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_bam_readcount_files_generates_expected_tsv(self):
        convert_vcf_input_file       = os.path.join(self.test_data_dir, 'full_input.vcf')
        convert_vcf_output_file      = tempfile.NamedTemporaryFile()
        convert_vcf_tdna_snvs_file   = os.path.join(self.test_data_dir, 'snvs.bam_readcount')
        convert_vcf_tdna_indels_file = os.path.join(self.test_data_dir, 'indels.bam_readcount')

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : convert_vcf_tdna_snvs_file,
            'tdna_indels_coverage_file'  : convert_vcf_tdna_indels_file,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_bam_readcount.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_multiple_transcripts_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_multiple_transcripts.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_multiple_transcripts.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_multiple_transcripts_per_alt_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_multiple_transcripts_per_alt.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_multiple_transcripts_per_alt.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_mutation_at_relative_beginning_of_full_sequence_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_mutation_at_relative_beginning_of_full_sequence.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_beginning_of_full_sequence.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_mutation_at_relative_end_of_full_sequence_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_mutation_at_relative_end_of_full_sequence.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_mutation_at_relative_end_of_full_sequence.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_position_out_of_bounds_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_position_out_of_bounds.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_position_out_of_bounds.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_short_wildtype_sequence_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_short_wildtype_sequence.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_short_wildtype_sequence.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_frameshift_variant_feature_elongation_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_elongation.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_elongation.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_frameshift_variant_feature_truncation_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_frameshift_variant_feature_truncation.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_frameshift_variant_feature_truncation.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_insertion_amino_acid_replacement_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_replacement.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_replacement.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_deletion_amino_acid_replacement_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_replacement.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_replacement.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_insertion_amino_acid_insertion_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_insertion_aa_insertion.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_insertion_aa_insertion.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_inframe_deletion_amino_acid_deletion_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_inframe_deletion_aa_deletion.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_inframe_deletion_aa_deletion.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))

    def test_input_vcf_with_uncalled_genotype_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_uncalled_genotype.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_uncalled_genotype.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))


    def test_input_vcf_with_hom_ref_genotype_generates_expected_tsv(self):
        convert_vcf_input_file  = os.path.join(self.test_data_dir, 'input_hom_ref_genotype.vcf')
        convert_vcf_output_file = tempfile.NamedTemporaryFile()

        convert_vcf_params = {
            'input_file'                 : convert_vcf_input_file,
            'output_file'                : convert_vcf_output_file.name,
            'gene_expn_file'             : None,
            'transcript_expn_file'       : None,
            'normal_snvs_coverage_file'  : None,
            'normal_indels_coverage_file': None,
            'tdna_snvs_coverage_file'    : None,
            'tdna_indels_coverage_file'  : None,
            'trna_snvs_coverage_file'    : None,
            'trna_indels_coverage_file'  : None,
        }
        converter = InputFileConverter(**convert_vcf_params)

        self.assertFalse(converter.execute())
        expected_output_file = os.path.join(self.test_data_dir, 'output_hom_ref_genotype.tsv')
        self.assertTrue(cmp(convert_vcf_output_file.name, expected_output_file))
