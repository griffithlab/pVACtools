import unittest
import os
import re
import sys
import tempfile
import csv
import py_compile
import shutil
import json
import requests
import time
from subprocess import run, PIPE, Popen, DEVNULL
from filecmp import cmp
pvac_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(pvac_dir)

#requests.get(, params={})
#requests.post

class APITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.urlBase = 'http://localhost:8080/api/v1'
        cls.pVac_directory =  pvac_dir
        cls.server_directory = os.path.join(
            pvac_dir,
            'pvacseq',
            'server'
        )
        cls.test_data_directory = os.path.join(
            cls.pVac_directory,
            'tests',
            'test_data',
            'api'
        )
        cls.api_process = None
        cls.restoreConfig = False
        cls.restoreData = False
        app_config_dir = os.path.expanduser(os.path.join(
            '~',
            '.pvacseq'
        ))
        app_data_dir = os.path.expanduser(os.path.join(
            '~',
            'pVAC-Seq'
        ))
        if os.path.isdir(app_config_dir):
            shutil.rmtree(app_config_dir+'.bak', ignore_errors=True)
            shutil.copytree(app_config_dir, app_config_dir+'.bak')
            shutil.rmtree(app_config_dir, ignore_errors=True)
            cls.restoreConfig = True
        if os.path.isdir(app_data_dir):
            shutil.rmtree(app_data_dir+'.bak', ignore_errors=True)
            shutil.copytree(app_data_dir, app_data_dir+'.bak')
            shutil.rmtree(app_data_dir, ignore_errors=True)
            cls.restoreData = True

    @classmethod
    def tearDownClass(cls):
        if isinstance(cls.api_process, Popen):
            cls.api_process.terminate()
        app_config_dir = os.path.expanduser(os.path.join(
            '~',
            '.pvacseq'
        ))
        app_data_dir = os.path.expanduser(os.path.join(
            '~',
            'pVAC-Seq'
        ))
        shutil.rmtree(app_config_dir, ignore_errors=True)
        shutil.rmtree(app_data_dir, ignore_errors=True)
        if cls.restoreConfig:
            shutil.copytree(app_config_dir+'.bak', app_config_dir)
            shutil.rmtree(app_config_dir+'.bak', ignore_errors=True)
        if cls.restoreData:
            shutil.copytree(app_data_dir+'.bak', app_data_dir)
            shutil.rmtree(app_data_dir+'.bak', ignore_errors=True)

    def test_app_compiles(self):
        compiled_app_path = py_compile.compile(os.path.join(
            self.server_directory,
            'app.py'
        ))
        self.assertTrue(compiled_app_path)
        controllers_dir = os.path.join(
            self.server_directory,
            'controllers'
        )
        for filename in os.listdir(controllers_dir):
            if filename.endswith('.py'):
                compiled_controller_path = py_compile.compile(os.path.join(
                    controllers_dir,
                    filename
                ))
                self.assertTrue(compiled_controller_path)

    def setUp(self):
        if APITests.api_process is None:
            print("Starting API...")
            APITests.api_process = Popen(
                [
                    'python',
                    '-m',
                    os.path.relpath(os.path.join(
                        self.server_directory,
                        'app'
                    ), pvac_dir).replace(os.sep, '.')
                ],
                # stdout = DEVNULL,
                # stderr = DEVNULL
            )
            time.sleep(5)

    def tearDown(self):
        time.sleep(2.5)

    def start_basic_run(self):
        response = requests.post(
            self.urlBase+'/staging',
            timeout = 5,
            data={
                'input':os.path.join(
                    self.test_data_directory,
                    'input.vcf'
                ),
                'samplename':'endpoint_input',
                'alleles':'HLA-G*01:09',
                'prediction_algorithms':'NetMHC',
            }
        )
        return response.json()

    def test_app_startup(self):
        #test config dir exists
        self.assertTrue(os.path.isdir(os.path.expanduser(os.path.join(
            '~',
            '.pvacseq'
        ))))
        #test data dirs
        self.assertTrue(os.path.isdir(os.path.expanduser(os.path.join(
            '~',
            'pVAC-Seq',
            'archive'
        ))))
        self.assertTrue(os.path.isdir(os.path.expanduser(os.path.join(
            '~',
            'pVAC-Seq',
            'input'
        ))))
        self.assertTrue(os.path.isdir(os.path.expanduser(os.path.join(
            '~',
            'pVAC-Seq',
            'results'
        ))))

    def test_endpoint_start(self):
        response = requests.post(
            self.urlBase + '/staging',
            timeout = 5,
            data={
                'input':os.path.join(
                    self.test_data_directory,
                    'input.vcf'
                ),
                'samplename':'endpoint_start',
                'alleles':'HLA-G*01:09,HLA-E*01:01',
                'prediction_algorithms':'NetMHC,PickPocket',
                'epitope_lengths':'9,10',
                'gene_expn_file':os.path.join(
                    self.test_data_directory,
                    'genes.fpkm_tracking'
                ),
                'transcript_expn_file':os.path.join(
                    self.test_data_directory,
                    'isoforms.fpkm_tracking'
                ),
                'tdna_snvs_coverage_file':os.path.join(
                    self.test_data_directory,
                    'snvs.bam_readcount'
                ),
                'tdna_indels_coverage_file':os.path.join(
                    self.test_data_directory,
                    'indels.bam_readcount'
                ),
                'top_score_metric':'lowest',
                'keep_tmp_files':'on',
                'netmhc_stab':'on',
                'net_chop_method':'cterm',
                'tdna_vaf':'40'

            }
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        self.assertTrue(re.match(r'\d+', response.content.decode()))

    def test_endpoint_input(self):
        shutil.copyfile(
            os.path.join(
                self.test_data_directory,
                'input.vcf'
            ),
            os.path.expanduser(os.path.join(
                '~',
                'pVAC-Seq',
                'input',
                'input.vcf'
            ))
        )
        response = requests.get(
            self.urlBase+'/input',
            timeout=5,
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        self.assertTrue(re.search(r'input\.vcf', response.content.decode()))
        input_manifest = response.json()
        vcf_id = list(filter(lambda x:x['name']=='input.vcf', input_manifest))[0]
        response = requests.post(
            self.urlBase+'/staging',
            timeout = 5,
            data={
                'input':vcf_id['fileID'],
                'samplename':'endpoint_input',
                'alleles':'HLA-G*01:09',
                'prediction_algorithms':'NetMHC',
            }
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        self.assertTrue(re.match(r'\d+', response.content.decode()))

    def test_endpoint_processes(self):
        response = requests.post(
            self.urlBase + '/staging',
            timeout = 5,
            data={
                'input':os.path.join(
                    self.test_data_directory,
                    'input.vcf'
                ),
                'samplename':'endpoint_processes',
                'alleles':'HLA-G*01:09',
                'prediction_algorithms':'NetMHC',
            }
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        self.assertTrue(re.match(r'\d+', response.content.decode()))
        processID = response.json()
        response = requests.get(
            self.urlBase+'/processes',
            timeout = 5,
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        targetResult = [item for item in response.json() if item['id'] == processID]
        self.assertTrue(len(targetResult))
        targetResult = targetResult[0]
        self.assertIsInstance(targetResult, dict)

        self.assertIn('attached', targetResult)
        self.assertIsInstance(targetResult['attached'], bool)

        self.assertIn('command', targetResult)
        self.assertIsInstance(targetResult['command'], str)
        self.assertTrue(targetResult['command'])

        self.assertIn('files', targetResult)
        self.assertIsInstance(targetResult['files'], list)
        self.assertTrue(targetResult['files'])
        self.assertIsInstance(targetResult['files'][0], dict)

        self.assertIn('id', targetResult)
        self.assertIsInstance(targetResult['id'], int)
        self.assertEqual(targetResult['id'], processID)

        self.assertIn('output', targetResult)
        self.assertIsInstance(targetResult['output'], str)
        self.assertTrue(targetResult['output'])

        self.assertIn('parameters', targetResult)
        self.assertIsInstance(targetResult['parameters'], dict)
        self.assertTrue(targetResult['parameters'])

        self.assertIn('pid', targetResult)
        self.assertIsInstance(targetResult['pid'], int)
        # self.assertGreater(targetResult['pid'], APITests.api_process.pid)

        self.assertIn('results_url', targetResult)
        self.assertIsInstance(targetResult['results_url'], str)
        self.assertTrue(targetResult['results_url'])

        self.assertIn('running', targetResult)
        self.assertIsInstance(targetResult['running'], bool)

        self.assertIn('url', targetResult)
        self.assertIsInstance(targetResult['url'], str)
        self.assertTrue(targetResult['url'])

    def test_endpoint_process_info(self):
        response = requests.get(
            self.urlBase + '/processes',
            timeout = 5
        )
        self.assertEqual(response.status_code,200)
        process_list = response.json()
        if not len(process_list):
            process_list = [{'id':self.start_basic_run()}]
        response = requests.get(
            self.urlBase + '/processes/%d'%process_list[0]['id'],
            timeout=5,
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        process_data = response.json()
        self.assertIsInstance(process_data, dict)

        self.assertIn('attached', process_data)
        self.assertIsInstance(process_data['attached'], bool)

        self.assertIn('command', process_data)
        self.assertIsInstance(process_data['command'], str)
        self.assertTrue(process_data['command'])

        self.assertIn('files', process_data)
        self.assertIsInstance(process_data['files'], list)
        self.assertTrue(process_data['files'])
        self.assertIsInstance(process_data['files'][0], dict)

        self.assertIn('id', process_data)
        self.assertEqual(process_data['id'], process_list[0]['id'])

        self.assertIn('log', process_data)
        self.assertIsInstance(process_data['log'], list)
        self.assertTrue(process_data['log'])

        self.assertIn('log_updated_at', process_data)
        self.assertIsInstance(process_data['log_updated_at'], int)
        self.assertLess(process_data['log_updated_at'], time.time())

        self.assertIn('output', process_data)
        self.assertIsInstance(process_data['output'], str)
        self.assertTrue(process_data['output'])

        self.assertIn('parameters', process_data)
        self.assertIsInstance(process_data['parameters'], dict)
        self.assertTrue(process_data['parameters'])

        self.assertIn('pid', process_data)
        self.assertIsInstance(process_data['pid'], int)
        # self.assertGreater(process_data['pid'], APITests.api_process.pid)

        self.assertIn('results_url', process_data)
        self.assertIsInstance(process_data['results_url'], str)
        self.assertTrue(process_data['results_url'])

        self.assertIn('running', process_data)
        self.assertIsInstance(process_data['running'], bool)

        self.assertIn('status', process_data)
        self.assertIsInstance(process_data['status'], str)
        self.assertTrue(process_data['status'])

    def test_endpoint_process_results(self):
        response = requests.get(
            self.urlBase + '/processes',
            timeout = 5
        )
        self.assertEqual(response.status_code,200)
        process_list = response.json()
        if not len(process_list):
            self.start_basic_run()
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
        process_list.sort(key = lambda x:x['id'])
        while process_list[-1]['running']:
            time.sleep(1)
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
            process_list.sort(key = lambda x:x['id'])
        response = requests.get(
            'http://localhost:8080'+process_list[-1]['results_url'],
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        results = response.json()
        self.assertIsInstance(results, list)
        for item in results:
            self.assertIsInstance(item, dict)

            self.assertIn('description', item)
            self.assertIsInstance(item['description'], str)
            self.assertTrue(item['description'])

            self.assertIn('display_name', item)
            self.assertIsInstance(item['display_name'], str)
            self.assertTrue(item['display_name'])

            self.assertIn('fileID', item)
            self.assertIsInstance(item['fileID'], str)
            self.assertGreaterEqual(int(item['fileID']), 0)

            self.assertIn('rows', item)
            self.assertIsInstance(item['rows'], int)
            self.assertGreaterEqual(item['rows'], -1)

            self.assertIn('size', item)
            self.assertIsInstance(item['size'], str)
            self.assertRegex(item['size'], r'\d+\.\d\d\d KB')

            self.assertIn('url', item)
            self.assertIsInstance(item['url'], str)
            self.assertTrue(item['url'])

    def test_endpoint_process_results_data(self):
        response = requests.get(
            self.urlBase + '/processes',
            timeout = 5
        )
        self.assertEqual(response.status_code,200)
        process_list = response.json()
        if not len(process_list):
            self.start_basic_run()
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
        process_list.sort(key = lambda x:x['id'])
        while process_list[-1]['running']:
            time.sleep(1)
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
            process_list.sort(key = lambda x:x['id'])
        response = requests.get(
            'http://localhost:8080'+process_list[-1]['results_url'],
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        results = response.json()
        for item in results:
            if item['display_name'].endswith('.tsv') and item['rows']>0:
                response = requests.get(
                    'http://localhost:8080'+item['url'],
                    timeout=5,
                )
                self.assertEqual(response.status_code, 200, response.content.decode())
                data = response.json()
                self.assertIsInstance(data, list)
                for row in data:
                    self.assertIsInstance(row, dict)
                    self.assertIn('rowid', row)
        #check the column name endpoint
        response = requests.get(
            'http://localhost:8080'+process_list[-1]['results_url']+'/cols',
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        #check the schema endpoint
        response = requests.get(
            'http://localhost:8080'+process_list[-1]['results_url']+'/schema',
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())

    def test_endpoint_stop(self):
        processID = self.start_basic_run()
        response = requests.get(
            self.urlBase+'/stop/%d'%processID,
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        response = requests.get(
            self.urlBase+'/processes/%d'%processID
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        data = response.json()
        self.assertIsInstance(data, dict)
        self.assertIn('running', data)
        self.assertFalse(data['running'])

    @unittest.skip("Shutdown is disabled")
    def test_endpoint_shutdown(self):
        processID = self.start_basic_run()
        response = requests.get(
            self.urlBase+'/shutdown',
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        data = response.json()
        self.assertIsInstance(data, list)
        self.assertIn(processID, data)

    def test_full_api(self):
        response = requests.post(
            self.urlBase + '/staging',
            timeout = 5,
            data={
                'input':os.path.join(
                    self.test_data_directory,
                    'input.vcf'
                ),
                'samplename':'endpoint_full',
                'alleles':'HLA-G*01:09,HLA-E*01:01',
                'prediction_algorithms':'NetMHC,PickPocket',
                'epitope_lengths':'9,10',
                'gene_expn_file':os.path.join(
                    self.test_data_directory,
                    'genes.fpkm_tracking'
                ),
                'transcript_expn_file':os.path.join(
                    self.test_data_directory,
                    'isoforms.fpkm_tracking'
                ),
                'tdna_snvs_coverage_file':os.path.join(
                    self.test_data_directory,
                    'snvs.bam_readcount'
                ),
                'tdna_indels_coverage_file':os.path.join(
                    self.test_data_directory,
                    'indels.bam_readcount'
                ),
                'top_score_metric':'lowest',
                'keep_tmp_files':'on',
                'netmhc_stab':'on',
                'net_chop_method':'cterm',
                'tdna_vaf':'40'

            }
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        processID = response.json()
        self.assertIsInstance(processID, int)
        response = requests.get(
            self.urlBase+'/processes/%d'%processID,
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        process_data = response.json()
        self.assertIsInstance(process_data, dict)
        self.assertIn('running', process_data)
        while process_data['running']:
            time.sleep(1)
            response = requests.get(
                self.urlBase+'/processes/%d'%processID,
                timeout = 5
            )
            self.assertEqual(response.status_code, 200, response.content.decode())
            process_data = response.json()
        self.assertIn('files', process_data)
        self.assertIsInstance(process_data['files'], list)
        finaltsv = [item for item in process_data['files'] if item['display_name'].endswith('.final.tsv')]
        self.assertTrue(finaltsv)
        finaltsv = finaltsv[0]
        raw_reader = open(os.path.join(
            self.test_data_directory,
            'Test.final.tsv'
        ))
        reader = csv.DictReader(raw_reader, delimiter='\t')
        response = requests.get(
            'http://localhost:8080'+finaltsv['url'],
            timeout=5,
            params={
                'count':10000
            }
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        content = response.json()
        self.assertIsInstance(content, list)
        output_file = tempfile.NamedTemporaryFile(mode='w+')
        writer = csv.DictWriter(
            output_file,
            fieldnames=reader.fieldnames,
            delimiter='\t',
            lineterminator='\n'
        )
        writer.writeheader()
        writer.writerows(content)
        output_file.flush()
        output_file.seek(0,0)
        raw_reader.seek(0,0)
        testlines = set(raw_reader.readlines())
        outputlines = set(raw_reader.readlines())
        self.assertEqual(len(testlines), len(outputlines), "Line count mismatch")
        self.assertFalse(testlines-outputlines, "Missing lines")
        self.assertFalse(outputlines-testlines, "Extra lines")
        self.assertFalse(set(raw_reader.readlines())^set(output_file.readlines()))
        raw_reader.close()
        output_file.close()

    def test_endpoint_allele(self):
        response = requests.get(
            self.urlBase+'/checkallele',
            timeout=5,
            params={
                'allele':'H2-IAb'
            }
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        self.assertTrue(response.json())
        response = requests.get(
            self.urlBase+'/checkallele',
            timeout=5,
            params={
                'allele':'NotARealAllele'
            }
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        self.assertFalse(response.json())

    @unittest.skip("Reset is disabled")
    def test_endpoing_reset(self):
        processID = self.start_basic_run()
        response = requests.get(
            self.urlBase+"/reset",
            timeout = 5,
            params={
                'clearall':1
            }
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        self.assertIn(processID, response.json())
        response = requests.get(
            self.urlBase+'/processes',
            timeout=5,
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        self.assertFalse(response.json())

    def test_endpoint_dropbox(self):
        shutil.copyfile(
            os.path.join(
                self.test_data_directory,
                'Test.final.tsv'
            ),
            os.path.expanduser(os.path.join(
                '~',
                'pVAC-Seq',
                'archive',
                'Test.final.tsv'
            ))
        )
        time.sleep(1)
        response = requests.get(
            self.urlBase+'/dropbox',
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        manifest = response.json()
        self.assertIsInstance(manifest, list)
        for item in manifest:
            self.assertIsInstance(item, dict)

            self.assertIn('description', item)
            self.assertIsInstance(item['description'], str)
            self.assertTrue(item['description'])

            self.assertIn('display_name', item)
            self.assertIsInstance(item['display_name'], str)
            self.assertTrue(item['display_name'])

            self.assertIn('fileID', item)
            self.assertIsInstance(item['fileID'], str)
            self.assertGreaterEqual(int(item['fileID']), 0)

            self.assertIn('rows', item)
            self.assertIsInstance(item['rows'], int)
            self.assertGreaterEqual(item['rows'], -1)

            self.assertIn('size', item)
            self.assertIsInstance(item['size'], str)
            self.assertRegex(item['size'], r'\d+\.\d\d\d KB')

            self.assertIn('url', item)
            self.assertIsInstance(item['url'], str)
            self.assertTrue(item['url'])

        target = [item for item in manifest if os.path.basename(item['display_name'])=='Test.final.tsv']
        self.assertTrue(target)
        target = target[0]
        self.assertEqual(target['rows'], 2)
        response = requests.get(
            'http://localhost:8080'+target['url'],
            timeout=5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        raw_reader = open(os.path.join(
            self.test_data_directory,
            'Test.final.tsv'
        ))
        reader = csv.DictReader(raw_reader, delimiter='\t')
        content = response.json()
        self.assertIsInstance(content, list)
        response = requests.get(
            'http://localhost:8080'+target['url']+'/cols',
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.content.decode())
        mapping = response.json()
        self.assertIsInstance(mapping, dict)
        output_file = tempfile.NamedTemporaryFile(mode='w+')
        writer = csv.DictWriter(
            output_file,
            fieldnames=reader.fieldnames,
            delimiter='\t',
            lineterminator='\n'
        )
        writer.writeheader()
        writer.writerows({mapping[key]:val for key, val in row.items() if key in mapping} for row in content)
        output_file.flush()
        output_file.seek(0,0)
        raw_reader.seek(0,0)
        testlines = set(raw_reader.readlines())
        outputlines = set(raw_reader.readlines())
        self.assertEqual(len(testlines), len(outputlines), "Line count mismatch")
        self.assertFalse(testlines-outputlines, "Missing lines")
        self.assertFalse(outputlines-testlines, "Extra lines")
        self.assertFalse(set(raw_reader.readlines())^set(output_file.readlines()))
        raw_reader.close()
        output_file.close()

    @unittest.skip("Filter is not implimented")
    def test_endpoint_filter(self):
        print("Heyo")
