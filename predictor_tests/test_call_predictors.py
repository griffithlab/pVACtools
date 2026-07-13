import unittest
import unittest.mock
import os
import sys
import re
import tempfile
from subprocess import call
from filecmp import cmp
import py_compile
import pandas as pd
from mock import patch

from pvactools.lib.call_predictors import CallPredictors
from pvactools.lib.prediction_class import PredictionClass, IEDB, NetMHCIIpan, NetMHCIIpanEL, NetMHCIIVersion

from tests.utils import *

def test_data_directory():
    return os.path.join(pvactools_directory(), 'predictor_tests', 'test_data')

class CallPredictorsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable_dir = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable     = os.path.join(cls.executable_dir, 'call_predictors.py')
        cls.test_data_dir  = test_data_directory()
        cls.additional_setup()

    @classmethod
    def additional_setup(cls):
        pass

class CallPredictorsCompileTests(CallPredictorsTests):
    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

class FilterResponseTests(CallPredictorsTests):
    def test_filter_response_ok(self):
        unfiltered_file = os.path.join(self.test_data_dir, 'unfiltered.txt')
        filtered_file   = os.path.join(self.test_data_dir, 'filtered.txt')
        with open(unfiltered_file, 'rb') as f:
            unfiltered_file_contents = f.read().rstrip()
        with open(filtered_file, 'rb') as f:
            filtered_file_contents = f.read().rstrip()
        filtered_response = IEDB.filter_response(unfiltered_file_contents)
        self.assertEqual(filtered_response, filtered_file_contents)
        filtered_response_on_filtered_file = IEDB.filter_response(filtered_file_contents)
        self.assertEqual(filtered_response_on_filtered_file, filtered_file_contents)

class CallClassIPredictorsTests(CallPredictorsTests):
    @classmethod
    def additional_setup(cls):
        cls.input_file     = os.path.join(cls.test_data_dir, 'input.fasta')
        cls.allele         = 'HLA-A*02:01'
        cls.epitope_length = 9
        cls.methods = ['ann', 'smmpmbec', 'smm', 'netmhcpan', 'netmhcpan_el']

    def test_iedb_methods_generate_expected_files(self):
        for method in self.methods:
            with tempfile.TemporaryDirectory() as output_dir:
                iedb_path = os.getenv('IEDB_PATH')
                if iedb_path is None:
                    raise Exception("IEDB_PATH env variable not set")

                class_name = PredictionClass.prediction_class_name_for_prediction_method(method)

                predictor_arguments = {
                    'input_file': self.input_file,
                    'sample_name': 'tmp',
                    'fasta_size': 800,
                    'allele': self.allele,
                    'epitope_length': self.epitope_length,
                    'prediction_algorithms': [class_name],
                    'iedb_executable_path': os.path.join(iedb_path, 'mhc_i', 'src', 'predict_binding.py'),
                    'iedb_retries': 5,
                    'n_threads': 1,
                    'output_dir': output_dir,
                }
                call_predictors = CallPredictors(**predictor_arguments)
                self.assertFalse(call_predictors.execute())

                expected_output_file = os.path.join(self.test_data_dir, 'output_%s.tsv' % method)
                self.assertTrue(cmp(call_predictors.output_files[0], expected_output_file))

    def test_mhcflurry_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': self.allele,
                'epitope_length': self.epitope_length,
                'prediction_algorithms': ['MHCflurry'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_mhcflurry.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[1,7,8])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[1,7,8])
            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_mhcnuggetsi_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': self.allele,
                'epitope_length': self.epitope_length,
                'prediction_algorithms': ['MHCnuggetsI'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_mhcnuggetsI.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,3,4,5])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[0,3,4,5])
            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_bigmhc_el_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': self.allele,
                'epitope_length': self.epitope_length,
                'prediction_algorithms': ['BigMHC_EL'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_bigmhc_el.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[1,5,6])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[1,5,6])
            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_bigmhc_im_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': self.allele,
                'epitope_length': self.epitope_length,
                'prediction_algorithms': ['BigMHC_IM'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_bigmhc_im.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[1,5,6])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[1,5,6])
            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_deepimmuno_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': self.allele,
                'epitope_length': self.epitope_length,
                'prediction_algorithms': ['DeepImmuno'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_deepimmuno.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,3,4])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[0,3,4])
            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_mixmhcpred_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': self.allele,
                'epitope_length': self.epitope_length,
                'prediction_algorithms': ['MixMHCpred'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_mixmhcpred.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,6,7])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[0,6,7])
            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_prime_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': self.allele,
                'epitope_length': self.epitope_length,
                'prediction_algorithms': ['PRIME'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_prime.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,8,9])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[0,8,9])
            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_prime_method_generates_expected_files_for_A02110(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': 'HLA-A*02:110',
                'epitope_length': self.epitope_length,
                'prediction_algorithms': ['PRIME'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_prime.A02110.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,8,9])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[0,8,9])
            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

class CallClassIIPredictorsTests(CallPredictorsTests):
    @classmethod
    def additional_setup(cls):
        cls.input_file     = os.path.join(cls.test_data_dir, 'input.15.fasta')
        cls.epitope_length = 15
        cls.allele         = 'H2-IAb'
        cls.methods = ['nn_align', 'netmhciipan_ba', 'netmhciipan_el', 'smm_align']

    def test_iedb_methods_generate_expected_files(self):
        for method in self.methods:
            with tempfile.TemporaryDirectory() as output_dir:
                iedb_path = os.getenv('IEDB_PATH')
                if iedb_path is None:
                    raise Exception("IEDB_PATH env variable not set")

                class_name = PredictionClass.prediction_class_name_for_prediction_method(method)

                predictor_arguments = {
                    'input_file': self.input_file,
                    'sample_name': 'tmp',
                    'fasta_size': 800,
                    'allele': self.allele,
                    'epitope_length': self.epitope_length,
                    'prediction_algorithms': [class_name],
                    'iedb_executable_path': os.path.join(iedb_path, 'mhc_ii', 'mhc_II_binding.py'),
                    'iedb_retries': 5,
                    'n_threads': 1,
                    'output_dir': output_dir,
                }
                call_predictors = CallPredictors(**predictor_arguments)
                self.assertFalse(call_predictors.execute())
                expected_output_file = os.path.join(self.test_data_dir, 'output_%s.tsv' % method)
                self.assertTrue(cmp(call_predictors.output_files[0], expected_output_file))

    def test_mhcnuggetsii_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': 'DPA1*01:03-DPB1*01:01',
                'epitope_length': self.epitope_length,
                'prediction_algorithms': ['MHCnuggetsII'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_mhcnuggetsII.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,2,3])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[0,2,3])
            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_netmhciipan_method_with_version(self):
        with tempfile.TemporaryDirectory() as output_dir:
            iedb_path = os.getenv('IEDB_PATH')
            if iedb_path is None:
                raise Exception("IEDB_PATH env variable not set")

            NetMHCIIVersion.netmhciipan_version = '4.2'

            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': 'DRB1*01:01',
                'epitope_length': 12,
                'prediction_algorithms': ['NetMHCIIpan'],
                'iedb_executable_path': os.path.join(iedb_path, 'mhc_ii', 'mhc_II_binding.py'),
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())

            expected_output_file = os.path.join(self.test_data_dir, 'output_netmhciipan-4.2.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,2,3])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[0,2,3])

            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_netmhciipan_el_method_with_version(self):
        with tempfile.TemporaryDirectory() as output_dir:
            iedb_path = os.getenv('IEDB_PATH')
            if iedb_path is None:
                raise Exception("IEDB_PATH env variable not set")

            NetMHCIIVersion.netmhciipan_version = '4.2'

            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': 'DRB1*01:01',
                'epitope_length': 12,
                'prediction_algorithms': ['NetMHCIIpanEL'],
                'iedb_executable_path': os.path.join(iedb_path, 'mhc_ii', 'mhc_II_binding.py'),
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())

            expected_output_file = os.path.join(self.test_data_dir, 'output_netmhciipan_el-4.2.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,2,3])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[0,2,3])

            pd.testing.assert_frame_equal(
                expected_df,
                actual_df,
                check_like=True,
                check_exact=False,
                rtol=1e-3,
                atol=2e-2
            )

    def test_mixmhc2pred_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': 'DRB1*04:05',
                'epitope_length': 12,
                'prediction_algorithms': ['MixMHC2pred'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_mixmhc2pred.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[0,13,14])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[0,13,14])
            pd.testing.assert_frame_equal(expected_df, actual_df, check_like=True, check_exact=False)

    def test_immuscope_im_method_generates_expected_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            predictor_arguments = {
                'input_file': self.input_file,
                'sample_name': 'tmp',
                'fasta_size': 800,
                'allele': 'DRB1*01:01',
                'epitope_length': 15,
                'prediction_algorithms': ['ImmuScope_IM'],
                'iedb_executable_path': None,
                'iedb_retries': 5,
                'n_threads': 1,
                'output_dir': output_dir,
            }
            call_predictors = CallPredictors(**predictor_arguments)
            self.assertFalse(call_predictors.execute())
            expected_output_file = os.path.join(self.test_data_dir, 'output_immuscope_im.tsv')
            expected_df = pd.read_csv(expected_output_file, sep="\t", index_col=[1,5,6])
            actual_df = pd.read_csv(call_predictors.output_files[0], sep="\t", index_col=[1,5,6])
            pd.testing.assert_frame_equal(expected_df, actual_df, check_like=True, check_exact=False)

if __name__ == '__main__':
    unittest.main()
