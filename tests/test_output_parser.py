import unittest
import os
import sys
import tempfile
import py_compile
import requests

from pvactools.lib.output_parser import PvacseqOutputParser, PvacbindOutputParser
from tests.utils import *

class OutputParserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        executable_dir    = os.path.join(pvactools_directory(), 'pvactools', 'lib')
        cls.executable    = os.path.join(executable_dir, 'output_parser.py')
        cls.test_data_dir = os.path.join(pvactools_directory(), 'tests', 'test_data', 'output_parser')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_parse_output_runs_and_produces_expected_output_with_all_class_i_files(self):
        prediction_files = [
            os.path.join(self.test_data_dir, "input.NetMHC.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.SMM.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.SMMPMBEC.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.NetMHCpan.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.NetMHCpanEL.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.BigMHC_EL.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.BigMHC_IM.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.DeepImmuno.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.MHCflurry.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.MHCnuggetsI.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.MixMHCpred.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.PRIME.HLA-A*02:01.9.1-800.tsv"),
        ]
        tsv_file = os.path.join(self.test_data_dir, "input.all_class_i.tsv")
        key_files = [os.path.join(self.test_data_dir, "input.all_class_i.1-800.key")]
        output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'prediction_files': prediction_files,
            'tsv_file'        : tsv_file,
            'key_files'       : key_files,
            'output_file'     : output_file.name,
            'sample_name'     : 'input',
            'flurry_state'    : 'both',
        }
        parser = PvacseqOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output.all_class_i.tsv")
        self.assertTrue(compare(output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_with_all_class_i_files_normalized_percentiles(self):
        prediction_files = [
            os.path.join(self.test_data_dir, "input.NetMHC.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.SMM.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.SMMPMBEC.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.NetMHCpan.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.NetMHCpanEL.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.BigMHC_EL.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.BigMHC_IM.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.DeepImmuno.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.MHCflurry.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.MHCnuggetsI.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.MixMHCpred.HLA-A*02:01.9.1-800.tsv"),
            os.path.join(self.test_data_dir, "input.PRIME.HLA-A*02:01.9.1-800.tsv"),
        ]
        tsv_file = os.path.join(self.test_data_dir, "input.all_class_i.tsv")
        key_files = [os.path.join(self.test_data_dir, "input.all_class_i.1-800.key")]
        output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'prediction_files'          : prediction_files,
            'tsv_file'                  : tsv_file,
            'key_files'                 : key_files,
            'output_file'               : output_file.name,
            'sample_name'               : 'input',
            'flurry_state'              : 'both',
            'use_normalized_percentiles': True,
        }
        parser = PvacseqOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output.all_class_i.normalized_percentiles.tsv")
        self.assertTrue(compare(output_file.name, expected_output_file))

    def test_parse_output_runs_and_produces_expected_output_with_all_class_ii_files(self):
        prediction_files = [
            os.path.join(self.test_data_dir, "HCC1395_TUMOR_DNA.MHCnuggetsII.DRB1*04:05.12.1-800.tsv"),
            os.path.join(self.test_data_dir, "HCC1395_TUMOR_DNA.MixMHC2pred.DRB1*04:05.12.1-800.tsv"),
            os.path.join(self.test_data_dir, "HCC1395_TUMOR_DNA.NetMHCIIpan.DRB1*04:05.12.1-800.tsv"),
            os.path.join(self.test_data_dir, "HCC1395_TUMOR_DNA.NetMHCIIpanEL.DRB1*04:05.12.1-800.tsv"),
            os.path.join(self.test_data_dir, "HCC1395_TUMOR_DNA.NNalign.DRB1*04:05.12.1-800.tsv"),
            os.path.join(self.test_data_dir, "HCC1395_TUMOR_DNA.SMMalign.DRB1*04:05.12.1-800.tsv"),
            os.path.join(self.test_data_dir, "HCC1395_TUMOR_DNA.ImmuScope_IM.DRB1*04:05.12.1-800.tsv"),
        ]
        tsv_file = os.path.join(self.test_data_dir, "input.all_class_ii.tsv")
        key_files = [os.path.join(self.test_data_dir, "input.all_class_ii.1-800.key")]
        output_file = tempfile.NamedTemporaryFile()

        parse_output_params = {
            'prediction_files': prediction_files,
            'tsv_file'        : tsv_file,
            'key_files'       : key_files,
            'output_file'     : output_file.name,
            'sample_name'     : 'input',
            'flurry_state'    : 'both',
        }
        parser = PvacseqOutputParser(**parse_output_params)

        self.assertFalse(parser.execute())
        expected_output_file  = os.path.join(self.test_data_dir, "output.all_class_ii.tsv")
        self.assertTrue(compare(output_file.name, expected_output_file))

    def test_get_scores_None_percentile(self):
        parse_output_params = {
            'prediction_files': [],
            'tsv_file'        : None,
            'key_files'       : [],
            'output_file'     : '',
            'sample_name'     : 'input_percentile_none',
            'flurry_state'    : None
        }
        parser = PvacseqOutputParser(**parse_output_params)

        line = {'allele': 'HLA-C*03:03', 'seq_num': '106', 'start': '7', 'end': '15', 'length': '9', 'peptide': 'FARGVAQPL', 'core': 'FARGVAQPL', 'icore': 'FARGVAQPL', 'ic50': '5.8', 'rank': 'None'}
        method = 'NetMHCpan'
        scores = parser.get_scores(line, method)
        expected_scores = {'NetMHCpan': {'ic50': 5.8, 'percentile': 'NA'}}
        self.assertEqual(scores, expected_scores)

    def test_get_scores_empty_percentile(self):
        parse_output_params = {
            'prediction_files': [],
            'tsv_file'        : None,
            'key_files'       : [],
            'output_file'     : '',
            'sample_name'     : 'input_percentile_none',
            'flurry_state'    : None
        }
        parser = PvacseqOutputParser(**parse_output_params)

        line = {'allele': 'HLA-C*15:05', 'peptide': 'QPKPVIDG', 'ic50': '28394.79812418208', 'percentile': '', 'mhcflurry_processing_score': '0.0385892167687416', 'mhcflurry_presentation_score': '0.0040710675126724', 'mhcflurry_presentation_percentile': '62.74467391304348', 'seq_num': '1', 'start': '1'}
        method = 'MHCflurry'
        scores = parser.get_scores(line, method)
        expected_scores = {'MHCflurry': {'ic50': 28394.79812418208, 'percentile': 'NA'}}
        self.assertEqual(scores, expected_scores)
