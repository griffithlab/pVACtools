import unittest
import tempfile
import py_compile
import shutil
from filecmp import cmp
import os
import sys
import re
from subprocess import PIPE
from subprocess import run as subprocess_run
import unittest.mock
from mock import patch
import argparse

from pvactools.tools.pvacvector import *
import pvactools.tools.pvacvector.main as pvacvector_main
from tests.utils import *

def make_response(data, path, test_name):
    filename = 'response_%s_%s_%s_%s.tsv' % (data['allele'], data['length'], data['method'], test_name)
    reader = open(os.path.join(
        path,
        filename
    ), mode='r')
    response_obj = lambda :None
    response_obj.status_code = 200
    response_obj.text = reader.read()
    reader.close()
    return response_obj

def make_predict_response(input_file, allele, length, path, test_name):
    file_parts = input_file.rsplit(os.sep, 6)
    clip_count = file_parts[1]
    spacer = file_parts[2]
    filename = f'response_{allele}_{length}_ann_{test_name}_{clip_count}_{spacer}.tsv'
    reader = open(os.path.join(
        path,
        filename
    ), mode='r')
    response = reader.read()
    reader.close()
    return (response, 'w')

def test_data_directory():
    base_dir = pvactools_directory()
    return os.path.join(base_dir, 'tests', 'test_data', 'pvacvector')

#python -m unittest tests/test_pvacvector.py
class TestPvacvector(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base_dir = pvactools_directory()
        cls.python = sys.executable
        cls.executable = os.path.join(cls.base_dir, 'pvactools', 'tools', 'pvacvector', 'run.py')
        cls.test_run_name = 'test_pvacvector_produces_expected_output'
        cls.test_data_dir = test_data_directory()
        cls.test_data_temp_dir = os.path.join(cls.test_data_dir, 'tmp')
        cls.input_tsv = os.path.join(cls.test_data_dir, 'input_parse_test_input.tsv')
        cls.input_vcf = os.path.join(cls.test_data_dir, 'input_parse_test_input.vcf')
        cls.input_file = os.path.join(cls.test_data_dir, 'Test.vector.results.input.fa')
        cls.method = 'NetMHC'
        cls.keep_tmp = 'True'
        cls.allele = 'HLA-A*02:01'
        cls.epitope_length = '8'

    def test_run_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_pvacvector_compiles(self):
        compiled_path = py_compile.compile(os.path.join(
            self.base_dir,
            'pvactools',
            'tools',
            'pvacvector',
            'main.py'
        ))
        self.assertTrue(compiled_path)

    def test_parser(self):
        parser = pvacvector_main.define_parser()
        self.assertEqual(type(parser), argparse.ArgumentParser)

    def test_pvacvector_commands(self):
        pvac_script_path = os.path.join(
            self.base_dir,
            'pvactools',
            'tools',
            'pvacvector',
            'main.py'
        )
        usage_search = re.compile(r"usage: ")
        for command in [
            "run",
            "visualize",
            "download_example_data",
            ]:
            result = subprocess_run([
                sys.executable,
                pvac_script_path,
                command,
                '-h'
            ], shell=False, stdout=PIPE)
            self.assertFalse(result.returncode)
            self.assertRegex(result.stdout.decode(), usage_search)

    def test_visualize_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.base_dir,
            'pvactools',
            "tools",
            "pvacvector",
            "visualize.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_visualize_runs(self):
        if 'DISPLAY' in os.environ.keys():
            input_file = os.path.join(self.test_data_dir, 'Test.vector.results.output.fa')
            output_dir = tempfile.TemporaryDirectory()
            visualize.main([input_file, output_dir.name])
            output_dir.cleanup()
        else:
            with self.assertRaises(Exception) as context:
                input_file = os.path.join(self.test_data_dir, 'Test.vector.results.output.fa')
                output_dir = tempfile.TemporaryDirectory()
                visualize.main([input_file, output_dir.name])
                output_dir.cleanup()

    def test_download_example_data_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.base_dir,
            'pvactools',
            "tools",
            "pvacvector",
            "download_example_data.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_download_example_data_runs(self):
        output_dir = tempfile.TemporaryDirectory()
        download_example_data.main([output_dir.name])
        output_dir.cleanup()

    def test_pvacvector_fa_input_runs_and_produces_expected_output(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data: make_response(
            data,
            test_data_directory(),
            'fa_input',
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory()

            run.main([
                self.input_file,
                self.test_run_name,
                self.allele,
                self.method,
                output_dir.name,
                '-e1', self.epitope_length,
                '--allow-n-peptide-exclusion', '0',
                '--percentile-threshold-strategy', 'exploratory',
                '--binding-percentile-threshold', '100',
                '-k',
                '--fasta-size', '600'
            ])

            #vaccine design algorithm producing correct output with fasta input
            self.assertTrue(cmp(
                os.path.join(output_dir.name, self.test_run_name + '_results.fa'),
                os.path.join(self.test_data_dir, "Test.vector.results.output.fa")
            ))
            self.assertTrue(cmp(
                os.path.join(output_dir.name, self.test_run_name + '_results.dna.fa'),
                os.path.join(self.test_data_dir, "Test.vector.results.output.dna.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, 'junctions.tsv'),
                os.path.join(self.test_data_dir, 'Test.vector.results.output.junctions.tsv')
            ))

            if 'DISPLAY' in os.environ.keys():
                image_out = os.path.join(output_dir.name, 'vector.png')
                #vaccine visualization producing image
                self.assertTrue(os.path.exists(image_out))
                self.assertTrue(os.stat(image_out).st_size > 0)

            output_dir.cleanup()

    def test_pvacvector_clipping(self):
        with patch('pvactools.lib.prediction_class.IEDB.predict', unittest.mock.Mock(side_effect = lambda input_file, allele, length, path, retries, tmp_dir=None, log_dir=None: make_predict_response(
            input_file,
            allele,
            length,
            test_data_directory(),
            'clipping',
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory()

            run.main([
                os.path.join(self.test_data_dir, 'input_parse_test_output.fa'),
                self.test_run_name,
                self.allele,
                self.method,
                output_dir.name,
                '-e1', self.epitope_length,
                '-k',
                '-b', '32000',
                '--max-clip-length', '2',
                '--allow-n-peptide-exclusion', '0',
                '--percentile-threshold-strategy', 'exploratory',
                '--binding-percentile-threshold', '100',
                '--spacers', 'None,AAY',
                '--fasta-size', '400'
            ])

            self.assertTrue(compare(
                os.path.join(output_dir.name, "0", "None", "junctions.tsv"),
                os.path.join(self.test_data_dir, "clipped.0.None.junctions.tsv")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "0", "AAY", "junctions.tsv"),
                os.path.join(self.test_data_dir, "clipped.0.AAY.junctions.tsv")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "1", "None", "junctions.tsv"),
                os.path.join(self.test_data_dir, "clipped.1.None.junctions.tsv")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "1", "AAY", "junctions.tsv"),
                os.path.join(self.test_data_dir, "clipped.1.AAY.junctions.tsv")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "2", "None", "junctions.tsv"),
                os.path.join(self.test_data_dir, "clipped.2.None.junctions.tsv")
            ))

            self.assertTrue(compare(
                os.path.join(output_dir.name, "0", "None", "tmp", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.0.None.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "0", "None", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.0.None.kmer.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "0", "AAY", "tmp", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.0.AAY.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "0", "AAY", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.0.AAY.kmer.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "1", "None", "tmp", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.1.None.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "1", "None", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.1.None.kmer.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "1", "AAY", "tmp", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.1.AAY.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "1", "AAY", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.1.AAY.kmer.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "2", "None", "tmp", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.2.None.fa")
            ))
            self.assertTrue(compare(
                os.path.join(output_dir.name, "2", "None", "test_pvacvector_produces_expected_output.8.fa"),
                os.path.join(self.test_data_dir, "clipped.2.None.kmer.fa")
            ))

            self.assertTrue(compare(
                os.path.join(output_dir.name, "test_pvacvector_produces_expected_output_results.fa"),
                os.path.join(self.test_data_dir, "clipped.result.fa")
            ))

            output_dir.cleanup()

    def test_pvacvector_percentile_threshold(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            test_data_directory(),
            'percentile_threshold',
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory()

            run.main([
                os.path.join(self.test_data_dir, 'input_parse_test_output.fa'),
                self.test_run_name,
                self.allele,
                self.method,
                output_dir.name,
                '-e1', self.epitope_length,
                '-k',
                '-b', '32000',
                '--binding-percentile-threshold', '80',
                '--max-clip-length', '0',
                '--allow-n-peptide-exclusion', '0',
                '--spacers', 'None',
            ])

            self.assertTrue(compare(
                os.path.join(output_dir.name, "0", "None", "junctions.tsv"),
                os.path.join(self.test_data_dir, "percentile_threshold.junctions.tsv")
            ))

            output_dir.cleanup()

    @unittest.skip("non-deterministic order difference in vector design")
    def test_pvacvector_remove_peptides(self):
        output_dir = tempfile.TemporaryDirectory()

        self.assertFalse(run.main([
            self.input_file,
            self.test_run_name,
            self.allele,
            self.method,
            output_dir.name,
            '-e1', self.epitope_length,
            '-k',
            '-b', '22000',
            '--percentile-threshold-strategy', 'exploratory',
            '--binding-percentile-threshold', '100',
            '--max-clip-length', '0',
            '--spacers', 'None',
        ]))

        #The result from this run can differ even with TEST_FLAG=1 set, causing these tests to fail
        #If we figure out how to make this run deterministic - reenable these comparison tests
        #self.assertTrue(compare(
        #    os.path.join(output_dir.name, "without_MT.CASP10.S654R", "test_pvacvector_produces_expected_output_results.fa"),
        #    os.path.join(self.test_data_dir, "without_MT.CASP10.S654R.test_pvacvector_produces_expected_output_results.fa")
        #))
        #self.assertTrue(compare(
        #    os.path.join(output_dir.name, "without_MT.FAT3.R4848T", "test_pvacvector_produces_expected_output_results.fa"),
        #    os.path.join(self.test_data_dir, "without_MT.FAT3.R4848T.test_pvacvector_produces_expected_output_results.fa")
        #))
        #self.assertTrue(compare(
        #    os.path.join(output_dir.name, "without_MT.PEX1.V356I", "test_pvacvector_produces_expected_output_results.fa"),
        #    os.path.join(self.test_data_dir, "without_MT.PEX1.V356I.test_pvacvector_produces_expected_output_results.fa")
        #))
        #self.assertTrue(compare(
        #    os.path.join(output_dir.name, "without_MT.POM121C.G3107R", "test_pvacvector_produces_expected_output_results.fa"),
        #    os.path.join(self.test_data_dir, "without_MT.POM121C.G3107R.test_pvacvector_produces_expected_output_results.fa")
        #))
        #self.assertTrue(compare(
        #    os.path.join(output_dir.name, "without_MT.PRDM15.G654W", "test_pvacvector_produces_expected_output_results.fa"),
        #    os.path.join(self.test_data_dir, "without_MT.PRDM15.G654W.test_pvacvector_produces_expected_output_results.fa")
        #))
        #self.assertTrue(compare(
        #    os.path.join(output_dir.name, "without_MT.SUMF2.G23A", "test_pvacvector_produces_expected_output_results.fa"),
        #    os.path.join(self.test_data_dir, "without_MT.SUMF2.G23A.test_pvacvector_produces_expected_output_results.fa")
        #))
        #self.assertTrue(compare(
        #    os.path.join(output_dir.name, "without_MT.TP53.R157H", "test_pvacvector_produces_expected_output_results.fa"),
        #    os.path.join(self.test_data_dir, "without_MT.TP53.R157H.test_pvacvector_produces_expected_output_results.fa")
        #))

        #self.assertFalse(os.path.exists(
        #    os.path.join(output_dir.name, "without_MT.ACSL3.S345N", "test_pvacvector_produces_expected_output_results.fa"),
        #))
        #self.assertFalse(os.path.exists(
        #    os.path.join(output_dir.name, "without_MT.DTX3L.G501R", "test_pvacvector_produces_expected_output_results.fa"),
        #))
        #self.assertFalse(os.path.exists(
        #    os.path.join(output_dir.name, "without_MT.NRCAM.P838H", "test_pvacvector_produces_expected_output_results.fa"),
        #))

        output_dir.cleanup()

    def test_prevent_clipping_best_peptide(self):
        with patch('pvactools.lib.prediction_class.IEDB.predict', unittest.mock.Mock(side_effect = lambda input_file, allele, length, path, retries, tmp_dir=None, log_dir=None: make_predict_response(
            input_file,
            allele,
            length,
            test_data_directory(),
            'prevent_clipping',
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory()
            input_file = os.path.join(self.test_data_dir, 'Test.vector.prevent_clipping_best_peptide.input.fa')

            with self.assertLogs(level='INFO') as log:
                run.main([
                    input_file,
                    self.test_run_name,
                    self.allele,
                    self.method,
                    output_dir.name,
                    '-e1', self.epitope_length,
                    '-b', '22000',
                    '--percentile-threshold-strategy', 'exploratory',
                    '--binding-percentile-threshold', '100',
                    '--spacers', 'None',
                    '--fasta-size', '1000',
                ])
                self.assertIn("INFO:root:Clipping 1 amino acids off the end of peptide MT.14.LGALS2.ENST00000215886.4.missense.132E/Q would clip the best peptide. Skipping.", log.output)
                self.assertIn("INFO:root:Clipping 2 amino acids off the start of peptide MT.20.PKDREJ.ENST00000253255.5.missense.1875T/I would clip the best peptide. Skipping.", log.output)
                self.assertIn("INFO:root:Clipping 2 amino acids off the end of peptide MT.14.LGALS2.ENST00000215886.4.missense.132E/Q would clip the best peptide. Skipping.", log.output)

                best_peptides = [
                    "LYYSYGLLHI",
                    "ARPPQQPVP",
                    "YQPCDDMDY",
                    "MVCELAGNL",
                    "NMSSFKLKQ",
                    "EMSHFEPNE",
                    "RSRTYDMDV",
                    "KTVTISCTG"
                ]
                with open(os.path.join(output_dir.name, "test_pvacvector_produces_expected_output_results.fa"), "r") as file:
                    file_content = file.read()
                    for best_peptide in best_peptides:
                        self.assertIn(best_peptide, file_content)

                output_dir.cleanup()
