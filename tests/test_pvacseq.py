import unittest
import unittest.mock
import os
import re
import sys
import tempfile
import py_compile
from subprocess import run, PIPE
from filecmp import cmp
import yaml
pvac_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(pvac_dir)
import pvacseq.lib
from pvacseq.lib.pipeline import *
import datetime

def compare(path1, path2):
    r1 = open(path1)
    r2 = open(path2)
    result = not len(set(r1.readlines())^set(r2.readlines()))
    r1.close()
    r2.close()
    return result

def make_response(data, files, path):
    if not files:
        if 'length' in data:
            filename = 'response_%s_%s_%s.tsv' % (data['allele'], data['length'], data['method'])
        else:
            filename = 'response_%s_%s.tsv' % (data['allele'], data['method'])
        reader = open(os.path.join(
            path,
            filename
        ), mode='r')
        response_obj = lambda :None
        response_obj.status_code = 200
        response_obj.text = reader.read()
        reader.close()
        return response_obj
    else:
        basefile = os.path.basename(data['configfile'])
        reader = open(os.path.join(
            path,
            'net_chop.html' if basefile == 'NetChop.cf' else 'Netmhcstab.html'
        ), mode='rb')
        response_obj = lambda :None
        response_obj.status_code = 200
        response_obj.content = reader.read()
        reader.close()
        return response_obj

def generate_class_i_call(method, allele, length, path, input_path):
    reader = open(os.path.join(
        input_path,
        "MHC_Class_I",
        "tmp",
        "Test_21.fa.split_1-48"
    ), mode='r')
    text = reader.read()
    reader.close()
    return unittest.mock.call('http://tools-cluster-interface.iedb.org/tools_api/mhci/', data={
        'sequence_text': ""+text,
        'method':        method,
        'allele':        allele,
        'length':        length,
    })

def generate_class_ii_call(method, allele, path, input_path):
    reader = open(os.path.join(
        input_path,
        "MHC_Class_II",
        "tmp",
        "Test_31.fa.split_1-48"
    ), mode='r')
    text = reader.read()
    reader.close()
    return unittest.mock.call('http://tools-cluster-interface.iedb.org/tools_api/mhcii/', data={
        'sequence_text': ""+text,
        'method':        method,
        'allele':        allele,
    })

class PVACTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pVac_directory =  pvac_dir
        cls.test_data_directory = os.path.join(
            cls.pVac_directory,
            'tests',
            'test_data',
            'pvacseq'
        )
        cls.methods = {
            'ann': {
                'HLA-E*01:01': [9, 10]
            },
            'pickpocket': {
                'HLA-G*01:09': [9, 10],
                'HLA-E*01:01': [9, 10],
            },
        }
        cls.request_mock = unittest.mock.Mock(side_effect = lambda url, data, files=None: make_response(
            data,
            files,
            cls.test_data_directory
        ))
        pvacseq.lib.call_iedb.requests.post = cls.request_mock

    def test_pvacseq_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            'pvacseq',
            "pvacseq.py"
        ))
        self.assertTrue(compiled_pvac_path)

    def test_pvacseq_commands(self):
        pvac_script_path = os.path.join(
            self.pVac_directory,
            'pvacseq',
            "pvacseq.py"
            )
        usage_search = re.compile(r"usage: ")
        for command in [
            "binding_filter",
            "coverage_filter",
            "run",
            "generate_protein_fasta",
            "install_vep_plugin",
            "download_example_data",
            "valid_alleles",
            "config_files",
            ]:
            result = run([
                sys.executable,
                pvac_script_path,
                command,
                '-h'
            ], shell=False, stdout=PIPE)
            self.assertFalse(result.returncode)
            self.assertRegex(result.stdout.decode(), usage_search)

    def test_main_compiles(self):
        compiled_main_path = py_compile.compile(os.path.join(
            self.pVac_directory,
            'pvacseq',
            "lib",
            "main.py"
        ))
        self.assertTrue(compiled_main_path)

    def test_pvacseq_pipeline(self):
        output_dir = tempfile.TemporaryDirectory(dir = self.test_data_directory)

        additional_input_files = tempfile.NamedTemporaryFile('w')
        additional_input_file_list = {
            'gene_expn_file': os.path.join(self.test_data_directory, 'genes.fpkm_tracking'),
            'transcript_expn_file': os.path.join(self.test_data_directory, 'isoforms.fpkm_tracking'),
            'tdna_snvs_coverage_file': os.path.join(self.test_data_directory, 'snvs.bam_readcount'),
            'tdna_indels_coverage_file': os.path.join(self.test_data_directory, 'indels.bam_readcount'),
        }
        yaml.dump(additional_input_file_list, additional_input_files, default_flow_style=False)

        pvacseq.lib.main.main([
            os.path.join(self.test_data_directory, "input.vcf"),
            'Test',
            'HLA-G*01:09,HLA-E*01:01',
            'NetMHC',
            'PickPocket',
            output_dir.name,
            '-e', '9,10',
            '-i', additional_input_files.name,
            '--top-score-metric=lowest',
            '--keep-tmp-files',
            '--net-chop-method', 'cterm',
            '--netmhc-stab',
            '--tdna-vaf', '20',
            '-d', 'full',
        ])
        pvacseq.lib.main.main([
            os.path.join(self.test_data_directory, "input.vcf"),
            'Test',
            'H2-IAb',
            'NNalign',
            output_dir.name,
            '-i', additional_input_files.name,
            '--top-score-metric=lowest',
            '--keep-tmp-files',
            '-d', 'full',
        ])

        for file_name in (
            'Test.tsv',
            'Test.tsv_1-24',
            'Test.combined.parsed.tsv',
            'Test.filtered.binding.tsv',
            'Test.filtered.coverage.tsv',
            'Test.chop.tsv',
            'Test.stab.tsv',
            'Test.final.tsv',
        ):
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', file_name)
            expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', file_name)
            self.assertTrue(cmp(output_file, expected_file))

        for file_name in (
            'Test_21.fa.split_1-48',
            'Test_21.fa.split_1-48.key',
        ):
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
            expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'tmp', file_name)
            self.assertTrue(cmp(output_file, expected_file))

        for file_name in (
            'Test.HLA-G*01:09.9.parsed.tsv_1-48',
            'Test.HLA-G*01:09.10.parsed.tsv_1-48',
            'Test.HLA-E*01:01.9.parsed.tsv_1-48',
            'Test.HLA-E*01:01.10.parsed.tsv_1-48',
        ):
            output_file   = os.path.join(output_dir.name, 'MHC_Class_I', 'tmp', file_name)
            expected_file = os.path.join(self.test_data_directory, 'MHC_Class_I', 'tmp', file_name)
            self.assertTrue(compare(output_file, expected_file))

        self.assertEqual(len(self.request_mock.mock_calls), 9)
        #Class I output files
        methods = self.methods
        for method in methods.keys():
            for allele in methods[method].keys():
                for length in methods[method][allele]:
                    self.request_mock.assert_has_calls([
                        generate_class_i_call(method, allele, length, self.test_data_directory, output_dir.name)
                    ])
                    output_file   = os.path.join(output_dir.name, "MHC_Class_I", "tmp", 'Test.%s.%s.%s.tsv_1-48' % (method, allele, length))
                    expected_file = os.path.join(self.test_data_directory, "MHC_Class_I", "tmp", 'Test.%s.%s.%s.tsv_1-48' % (method, allele, length))
                    self.assertTrue(cmp(output_file, expected_file, False))

        #Class II output files
        for file_name in (
            'Test.tsv',
            'Test.tsv_1-24',
            'Test.combined.parsed.tsv',
            'Test.filtered.binding.tsv',
            'Test.filtered.coverage.tsv',
            'Test.final.tsv',
        ):
            output_file   = os.path.join(output_dir.name, 'MHC_Class_II', file_name)
            expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', file_name)
            self.assertTrue(cmp(output_file, expected_file, False))

        for file_name in (
            'Test_31.fa.split_1-48',
            'Test_31.fa.split_1-48.key',
            'Test.nn_align.H2-IAb.tsv_1-48',
        ):
            output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'tmp', file_name)
            expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', 'tmp', file_name)
            self.assertTrue(cmp(output_file, expected_file, False))

        for file_name in (
            'Test.H2-IAb.parsed.tsv_1-48',
        ):
            output_file   = os.path.join(output_dir.name, 'MHC_Class_II', 'tmp', file_name)
            expected_file = os.path.join(self.test_data_directory, 'MHC_Class_II', 'tmp', file_name)
            self.assertTrue(compare(output_file, expected_file))

        self.request_mock.assert_has_calls([
            generate_class_ii_call('nn_align', 'H2-IAb', self.test_data_directory, output_dir.name)
        ])

        output_dir.cleanup()

    def test_pvacseq_pipeline_sleep(self):
        output_dir_1 = tempfile.TemporaryDirectory()
        params_1 = [
            os.path.join(self.test_data_directory, "input.vcf"),
            'Test',
            'HLA-E*01:01',
            'NetMHC',
            output_dir_1.name,
            '-e', '9,10',
        ]
        os.environ["TEST_FLAG"] = '0'
        start_1 = datetime.datetime.now()
        pvacseq.lib.main.main(params_1)
        end_1 = datetime.datetime.now()
        duration_1 = (end_1 - start_1).total_seconds()
        output_dir_1.cleanup()

        output_dir_2 = tempfile.TemporaryDirectory()
        params_2 = [
            os.path.join(self.test_data_directory, "input.vcf"),
            'Test',
            'HLA-E*01:01',
            'NetMHC',
            output_dir_2.name,
            '-e', '9,10',
        ]
        os.environ["TEST_FLAG"] = '1'
        start_2 = datetime.datetime.now()
        pvacseq.lib.main.main(params_2)
        end_2 = datetime.datetime.now()
        duration_2 = (end_2 - start_2).total_seconds()
        output_dir_2.cleanup()

        self.assertTrue(duration_1 > duration_2)

