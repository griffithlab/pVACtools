import unittest
import os
import tempfile
from subprocess import call
from filecmp import cmp

class GenerateVariantSequences(unittest.TestCase):
    def runTest(self):
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        executable_dir = os.path.join(base_dir, 'pvac_seq')
        test_data_dir  = os.path.join(base_dir, 'test_data')

        sample_name             = 'Test'
        peptide_sequence_length = 21

        generate_variant_sequences_input_file  = os.path.join(test_data_dir, 'annotated_variants.tsv')
        generate_variant_sequences_output_file = tempfile.NamedTemporaryFile().name
        generate_variant_sequences_executable  = os.path.join(executable_dir, 'generate_variant_sequences.py')

        generate_variant_sequences_command = "python %s %s %s %s" % (generate_variant_sequences_executable, generate_variant_sequences_input_file, peptide_sequence_length, generate_variant_sequences_output_file)

        self.assertFalse(call(generate_variant_sequences_command, shell=True), 'GenerateVariantSequences command executes successfully')
        expected_output_file = os.path.join(test_data_dir, ("%s_%s.fa" % (sample_name, peptide_sequence_length)))
        self.assertTrue(cmp(generate_variant_sequences_output_file, expected_output_file), 'GenerateVariantSequences output as expected')

if __name__ == '__main__':
    unittest.main()
