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
        cls.input_n_mer = '25'

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
            "valid_alleles",
            "valid_algorithms",
            "allele_specific_cutoffs",
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

    def test_allele_specific_cutoffs_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.base_dir,
            'pvactools',
            "tools",
            "pvacvector",
            "allele_specific_cutoffs.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_allele_specific_cutoffs_runs(self):
        allele_specific_cutoffs.main([])

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

    def test_valid_alleles_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.base_dir,
            'pvactools',
            "tools",
            "pvacvector",
            "valid_alleles.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_valid_alleles_runs(self):
        valid_alleles.main(["-p", "SMM"])

    def test_valid_algorithms_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            self.base_dir,
            'pvactools',
            "tools",
            "pvacvector",
            "valid_algorithms.py"
        ))
        self.assertTrue(compiled_run_path)

    def test_valid_algorithms_runs(self):
        valid_algorithms.main("")

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
                '-n', self.input_n_mer,
                '--allow-n-peptide-exclusion', '0',
                '-k'
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

    def test_pvacvector_generate_fa_runs_and_produces_expected_output(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            test_data_directory(),
            'generate_fa',
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory()

            run.main([
                self.input_tsv,
                self.test_run_name,
                self.allele,
                self.method,
                output_dir.name,
                '-v', self.input_vcf,
                '-e1', self.epitope_length,
                '-n', self.input_n_mer,
                '--allow-n-peptide-exclusion', '0',
                '-k',
            ])

            #conversion from vcf to fasta file producing correct output, input file for vaccine design algorithm
            self.assertTrue(compare(
                    os.path.join(output_dir.name, "vector_input.fa"),
                    os.path.join(self.test_data_dir, "input_parse_test_output.fa")
                    ))


            if 'DISPLAY' in os.environ.keys():
                image_out = os.path.join(output_dir.name, 'vector.png')
                #vaccine visualization producing image
                self.assertTrue(os.path.exists(image_out))
                self.assertTrue(os.stat(image_out).st_size > 0)

            output_dir.cleanup()

    def test_pvacvector_generate_fa_with_epitope_at_beginning_of_transcript(self):
        with patch('requests.post', unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            test_data_directory(),
            'negative_start',
        ))) as mock_request:
            output_dir = tempfile.TemporaryDirectory()

            run.main([
                os.path.join(self.test_data_dir, "input_negative_start.tsv"),
                'H_MT-10109-005',
                self.allele,
                self.method,
                output_dir.name,
                '-v', os.path.join(self.test_data_dir, "input_negative_start.vcf.gz"),
                '-e1', self.epitope_length,
                '-n', self.input_n_mer,
                '--allow-n-peptide-exclusion', '0',
                '-k',
            ])

            self.assertTrue(compare(
                os.path.join(output_dir.name, "vector_input.fa"),
                os.path.join(self.test_data_dir, "output_negative_start.fa")
            ))
            output_dir.cleanup()

    def test_pvacvector_clipping(self):
        output_dir = tempfile.TemporaryDirectory()

        run.main([
            self.input_tsv,
            self.test_run_name,
            self.allele,
            self.method,
            output_dir.name,
            '-v', self.input_vcf,
            '-e1', self.epitope_length,
            '-n', self.input_n_mer,
            '-k',
            '-b', '32000',
            '--max-clip-length', '2',
            '--allow-n-peptide-exclusion', '0',
            '--spacers', 'None,AAY',
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
            os.path.join(output_dir.name, "0", "None", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv"),
            os.path.join(self.test_data_dir, "clipped.0.None.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "0", "None", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv.key"),
            os.path.join(self.test_data_dir, "clipped.0.None.fa.key")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "0", "AAY", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv"),
            os.path.join(self.test_data_dir, "clipped.0.AAY.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "0", "AAY", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv.key"),
            os.path.join(self.test_data_dir, "clipped.0.AAY.fa.key")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "1", "None", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv"),
            os.path.join(self.test_data_dir, "clipped.1.None.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "1", "None", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv.key"),
            os.path.join(self.test_data_dir, "clipped.1.None.fa.key")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "1", "AAY", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv"),
            os.path.join(self.test_data_dir, "clipped.1.AAY.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "1", "AAY", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv.key"),
            os.path.join(self.test_data_dir, "clipped.1.AAY.fa.key")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "2", "None", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv"),
            os.path.join(self.test_data_dir, "clipped.2.None.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "2", "None", "MHC_Class_I", "tmp", "test_pvacvector_produces_expected_output.fa.split_1-2.8.tsv.key"),
            os.path.join(self.test_data_dir, "clipped.2.None.fa.key")
        ))

        self.assertTrue(compare(
            os.path.join(output_dir.name, "test_pvacvector_produces_expected_output_results.fa"),
            os.path.join(self.test_data_dir, "clipped.result.fa")
        ))

        output_dir.cleanup()

    def test_pvacvector_percentile_threshold(self):
        output_dir = tempfile.TemporaryDirectory()

        run.main([
            self.input_tsv,
            self.test_run_name,
            self.allele,
            self.method,
            output_dir.name,
            '-v', self.input_vcf,
            '-e1', self.epitope_length,
            '-n', self.input_n_mer,
            '-k',
            '-b', '32000',
            '--percentile-threshold', '80',
            '--max-clip-length', '0',
            '--allow-n-peptide-exclusion', '0',
            '--spacers', 'None',
        ])

        self.assertTrue(compare(
            os.path.join(output_dir.name, "0", "None", "junctions.tsv"),
            os.path.join(self.test_data_dir, "percentile_threshold.junctions.tsv")
        ))

        output_dir.cleanup()

    def test_pvacvector_remove_peptides(self):
        output_dir = tempfile.TemporaryDirectory()

        run.main([
            self.input_file,
            self.test_run_name,
            self.allele,
            self.method,
            output_dir.name,
            '-e1', self.epitope_length,
            '-n', self.input_n_mer,
            '-k',
            '-b', '22000',
            '--max-clip-length', '0',
            '--spacers', 'None',
        ])

        self.assertTrue(compare(
            os.path.join(output_dir.name, "without_MT.CASP10.S654R", "test_pvacvector_produces_expected_output_results.fa"),
            os.path.join(self.test_data_dir, "without_MT.CASP10.S654R.test_pvacvector_produces_expected_output_results.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "without_MT.FAT3.R4848T", "test_pvacvector_produces_expected_output_results.fa"),
            os.path.join(self.test_data_dir, "without_MT.FAT3.R4848T.test_pvacvector_produces_expected_output_results.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "without_MT.PEX1.V356I", "test_pvacvector_produces_expected_output_results.fa"),
            os.path.join(self.test_data_dir, "without_MT.PEX1.V356I.test_pvacvector_produces_expected_output_results.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "without_MT.POM121C.G3107R", "test_pvacvector_produces_expected_output_results.fa"),
            os.path.join(self.test_data_dir, "without_MT.POM121C.G3107R.test_pvacvector_produces_expected_output_results.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "without_MT.PRDM15.G654W", "test_pvacvector_produces_expected_output_results.fa"),
            os.path.join(self.test_data_dir, "without_MT.PRDM15.G654W.test_pvacvector_produces_expected_output_results.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "without_MT.SUMF2.G23A", "test_pvacvector_produces_expected_output_results.fa"),
            os.path.join(self.test_data_dir, "without_MT.SUMF2.G23A.test_pvacvector_produces_expected_output_results.fa")
        ))
        self.assertTrue(compare(
            os.path.join(output_dir.name, "without_MT.TP53.R157H", "test_pvacvector_produces_expected_output_results.fa"),
            os.path.join(self.test_data_dir, "without_MT.TP53.R157H.test_pvacvector_produces_expected_output_results.fa")
        ))

        self.assertFalse(os.path.exists(
            os.path.join(output_dir.name, "without_MT.ACSL3.S345N", "test_pvacvector_produces_expected_output_results.fa"),
        ))
        self.assertFalse(os.path.exists(
            os.path.join(output_dir.name, "without_MT.DTX3L.G501R", "test_pvacvector_produces_expected_output_results.fa"),
        ))
        self.assertFalse(os.path.exists(
            os.path.join(output_dir.name, "without_MT.NRCAM.P838H", "test_pvacvector_produces_expected_output_results.fa"),
        ))

        output_dir.cleanup()

    def test_prevent_clipping_best_peptide(self):
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
                '-n', self.input_n_mer,
                '-k',
                '-b', '22000',
                '--spacers', 'None',
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
