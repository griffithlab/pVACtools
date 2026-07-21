import unittest
import unittest.mock
import os
import re
import sys
import py_compile
from subprocess import PIPE
from subprocess import run as subprocess_run
from filecmp import cmp
from mock import patch
import argparse
import logging
from testfixtures import LogCapture, StringComparison as S

from pvactools.lib.fasta_to_kmers import SequenceFastaToKmers
import pvactools.tools.pvacbind.main as pvacbind_main
from pvactools.tools.pvacbind import run
from tests.utils import *

def test_data_directory():
    return os.path.join(
        pvactools_directory(),
        'tests',
        'test_data',
        'pvacbind'
    )

class PvacbindTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pvactools_directory = pvactools_directory()
        cls.test_data_directory = test_data_directory()
        cls.methods = {
            'ann': {
                'HLA-E*01:01': [9, 10]
            },
            'pickpocket': {
                'HLA-G*01:09': [9, 10],
                'HLA-E*01:01': [9, 10],
            },
        }
        cls.peptide_fasta = os.path.join(pvactools_directory(), "tests", "test_data", "Homo_sapiens.GRCh38.pep.short.fa.gz")

    def test_pvacbind_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            'tools',
            'pvacbind',
            "main.py"
        ))
        self.assertTrue(compiled_pvac_path)

    def test_parser(self):
        parser = pvacbind_main.define_parser()
        self.assertEqual(type(parser), argparse.ArgumentParser)

    def test_pvacbind_commands(self):
        pvac_script_path = os.path.join(
            self.pvactools_directory,
            "pvactools",
            'tools',
            'pvacbind',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "run",
            ]:
            result = subprocess_run([
                sys.executable,
                pvac_script_path,
                command,
                '-h'
            ], shell=False, stdout=PIPE)
            self.assertFalse(result.returncode)
            self.assertRegex(result.stdout.decode(), usage_search)

    def test_run_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.pvactools_directory,
            "pvactools",
            "tools",
            "pvacbind",
            "run.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_process_stops(self):
        output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
        fasta_to_kmer_arguments = {
            'fasta': os.path.join(self.test_data_directory, "input_with_stops.fasta"),
            'output_dir': output_dir.name,
            'epitope_length': 9,
            'sample_name': 'Test',
        }
        SequenceFastaToKmers(**fasta_to_kmer_arguments).execute()
        output_file   = os.path.join(output_dir.name, 'Test.9.fa')
        expected_file = os.path.join(self.test_data_directory, 'output_with_stops.fasta')
        self.assertTrue(cmp(output_file, expected_file))
        output_dir.cleanup()

    def test_pvacbind_pipeline(self):
        with patch('pvactools.lib.call_predictors.requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            test_data_directory()
        ))) as mock_request, patch('pvactools.lib.net_chop.NetChop.post_query', unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            'net_chop.html'
        ))), patch('pvactools.lib.netmhc_stab.NetMHCStab.query_netmhcstabpan_server',  unittest.mock.Mock(side_effect = lambda url, data, timeout, files=None: mock_netchop_netmhcstabpan(
            data,
            files,
            self.test_data_directory,
            'Netmhcstab.html'
        ))):
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

            run.main([
                os.path.join(self.test_data_directory, "input.fasta"),
                'sample.name',
                'HLA-G*01:09,HLA-E*01:01,DRB1*11:01',
                'NetMHC',
                'PickPocket',
                'NNalign',
                output_dir.name,
                '-e1', '9,10',
                '-e2', '15',
                '--top-score-metric=lowest',
                '--top-score-metric2=ic50',
                '--keep-tmp-files',
                '--net-chop-method', 'cterm',
                '--netmhc-stab',
                '--run-reference-proteome-similarity',
                '--peptide-fasta', self.peptide_fasta,
                '--fasta-size', '3000',
            ])

            close_mock_fhs()

            #Shared output files
            for file_name in (
                'sample.name.9.fa',
                'sample.name.10.fa',
                'sample.name.15.fa',
            ):
                output_file   = os.path.join(output_dir.name, file_name)
                expected_file = os.path.join(self.test_data_directory, "run", file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'inputs.yml',
            ):
                output_file   = os.path.join(output_dir.name, 'log', file_name)
                self.assertTrue(os.path.exists(output_file))

            #Class I output files
            for file_name in (
                'sample.name.MHC_I.all_epitopes.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'run', 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.MHC_I.filtered.tsv',
                'sample.name.MHC_I.all_epitopes.aggregated.tsv',
                'sample.name.MHC_I.all_epitopes.aggregated.tsv.reference_matches',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
                expected_file = os.path.join(self.test_data_directory, 'run', 'MHC_Class_I', file_name.replace('sample.name', 'Test'))
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.HLA-G*01:09.9.parsed.tsv',
                'sample.name.HLA-E*01:01.9.parsed.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', '9', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'run', 'MHC_Class_I', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.HLA-G*01:09.10.parsed.tsv',
                'sample.name.HLA-E*01:01.10.parsed.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_I', '10', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'run', 'MHC_Class_I', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            #Class II output files
            for file_name in (
                'sample.name.MHC_II.all_epitopes.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', file_name)
                expected_file = os.path.join(self.test_data_directory, 'run', 'MHC_Class_II', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.MHC_II.filtered.tsv',
                'sample.name.MHC_II.all_epitopes.aggregated.tsv',
                'sample.name.MHC_II.all_epitopes.aggregated.tsv.reference_matches',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', file_name)
                expected_file = os.path.join(self.test_data_directory, 'run', 'MHC_Class_II', file_name.replace('sample.name', 'Test'))
                self.assertTrue(cmp(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            for file_name in (
                'sample.name.DRB1*11:01.15.parsed.tsv',
            ):
                output_file   = os.path.join(output_dir.name, 'MHC_Class_II', '15', 'tmp', file_name)
                expected_file = os.path.join(self.test_data_directory, 'run', 'MHC_Class_II', 'tmp', file_name.replace('sample.name', 'Test'))
                self.assertTrue(compare(output_file, expected_file), "files don't match %s - %s" %(output_file, expected_file))

            with self.assertRaises(SystemExit) as cm:
                run.main([
                    os.path.join(self.test_data_directory, "input.fasta"),
                    'sample.name',
                    'HLA-G*01:09,HLA-E*01:01,DRB1*11:01',
                    'NetMHC',
                    'PickPocket',
                    'NNalign',
                    output_dir.name,
                    '-e1', '9,10',
                    '-e2', '15',
                    '--top-score-metric2=ic50',
                    '--keep-tmp-files',
                    '--net-chop-method', 'cterm',
                    '--netmhc-stab',
                    '--run-reference-proteome-similarity',
                    '--peptide-fasta', self.peptide_fasta,
                    '--fasta-size', '3000',
                ])
            self.assertEqual(
                str(cm.exception),
                "Restart inputs are different from past inputs: \n" +
                "Past input: top_score_metric - lowest\n" +
                "Current input: top_score_metric - median\nAborting."
            )

            output_dir.cleanup()

    def test_duplicate_fasta_header(self):
        with self.assertRaises(Exception) as cm:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
            run.main([
                os.path.join(self.test_data_directory, "input.duplicate_header.fasta"),
                'Test',
                'HLA-A*02:01',
                'NetMHC',
                output_dir.name,
                '-e1', '8'
            ])
        self.assertEqual(
            str(cm.exception),
            'Duplicate fasta header "1". Please ensure that the input FASTA uses unique headers.'
        )
        output_dir.cleanup()

    def test_unsupported_amino_acid(self):
        logging.disable(logging.NOTSET)
        with LogCapture() as l:
            output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)
            fasta_to_kmer_arguments = {
                'fasta': os.path.join(self.test_data_directory, "input.unsupported_amino_acid.fasta"),
                'output_dir': output_dir.name,
                'epitope_length': 8,
                'sample_name': 'Test',
            }
            SequenceFastaToKmers(**fasta_to_kmer_arguments).execute()
            l.check_present(('root', 'WARNING', S("Record LPZLPPPP contains unsupported amino acids. Skipping.")))
            output_dir.cleanup()
