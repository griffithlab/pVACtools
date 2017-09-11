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
import signal
import time
import random
import postgresql as psql
from socket import *
from . import mock_api
from subprocess import run, PIPE, Popen, DEVNULL, TimeoutExpired
from filecmp import cmp
pvac_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(pvac_dir)

def parsedata(data):
    if re.match(r'^\d+\.\d+$', data):
        return float(data)
    if re.match(r'^\d+$', data):
        return int(data)
    if data == 'NA':
        return None
    return data

# checks if server has been successfully started
def check_opened():
    s = socket(AF_INET, SOCK_STREAM, 0)
    result = s.connect_ex(("127.0.0.1", 8080))
    if result == 0:
        s.close()
        return True
    s.close()
    return False

class APITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed()
        cls.urlBase = 'http://localhost:8080/api/v1'
        cls.pVac_directory =  pvac_dir
        cls.server_directory = os.path.join(
            pvac_dir,
            'utils',
            'pvacapi'
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
            os.kill(cls.api_process.pid, signal.SIGINT)
            try:
                cls.api_process.wait(1.5)
            except TimeoutExpired:
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
        db = psql.open("localhost/pvacseq")
        for row in db.prepare("SELECT table_name FROM information_schema.tables WHERE table_name LIKE 'data\__%\__%'")():
            name = row[0]
            if re.match(r'data_(dropbox|\d+)_\d+', name):
                print("DROP TABLE", name)
                db.execute("DROP TABLE %s"%name)

    def setUp(self):
        if APITests.api_process is None:
            print("Starting API...")
            APITests.api_process = Popen(
                [
                    sys.executable,
                    os.path.join(os.path.dirname(__file__),'mock_api.py'),
                    'api'
                ],
                stdout = DEVNULL,
                stderr = DEVNULL,
            )
            time.sleep(5)
            total_time = 5
            while not check_opened():
                if (total_time > 30):
                    break
                time.sleep(5)
                total_time += 5
            requests.get(
                self.urlBase+'/processes',
                timeout = 10
            )

    def tearDown(self):
        time.sleep(.5)

    def start_basic_run(self):
        response = requests.post(
            self.urlBase+'/staging',
            timeout = 5,
            data={
                'input':os.path.join(
                    self.test_data_directory,
                    'input.vcf'
                ),
                'samplename':'basic_run',
                'alleles':'HLA-E*01:01',
                'prediction_algorithms':'NetMHC',
                'force':True
            }
        )
        self.assertEqual(response.status_code, 201)
        return response.json()

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
            'dropbox'
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
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        self.assertTrue(re.search(r'input\.vcf', response.content.decode()))
        input_manifest = response.json()
        vcf_id = list(filter(lambda x:x['display_name']=='input.vcf', input_manifest))[0]
        response = requests.post(
            self.urlBase+'/staging',
            timeout = 5,
            data={
                'input':vcf_id['fileID'],
                'samplename':'endpoint_input',
                'alleles':'HLA-G*01:09',
                'prediction_algorithms':'NetMHC',
                'force':True
            }
        )
        self.assertEqual(response.status_code, 201, response.url+' : '+response.content.decode())
        result = response.json()
        self.assertTrue(re.match(r'\d+', str(result['processid'])))
        self.assertTrue(re.match(r'\d+', str(result['code'])))
        self.assertTrue(re.match(r'\S+', result['message']))
    
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
                'force':True
            }
        )
        self.assertEqual(response.status_code, 201, response.url+' : '+response.content.decode())
        result = response.json()
        self.assertTrue(re.match(r'\d+', str(result['processid'])))
        self.assertTrue(re.match(r'\d+', str(result['code'])))
        self.assertTrue(re.match(r'\S+', result['message']))     
        processID = result['processid']
        response = requests.get(
            self.urlBase+'/processes',
            timeout = 5,
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        targetResult = [item for item in response.json()['result'] if item['id'] == processID]
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

        self.assertIn('returncode', targetResult)
        self.assertIsInstance(targetResult['returncode'], int)

        self.assertIn('status', targetResult)
        self.assertIsInstance(targetResult['status'], int)

    def test_endpoint_process_info(self):
        response = requests.get(
            self.urlBase + '/processes',
            timeout = 5
        )
        self.assertEqual(response.status_code,200)
        process_list = response.json()
        if not len(process_list['result']):
            process_list = [{'id':self.start_basic_run()}]
        response = requests.get(
            self.urlBase + '/processes/%d'%process_list['result'][0]['id'],
            timeout=5,
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
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
        self.assertEqual(process_data['id'], process_list['result'][0]['id'])

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

        self.assertIn('last_message', process_data)
        self.assertIsInstance(process_data['last_message'], str)
        self.assertTrue(process_data['last_message'])

        self.assertIn('running', process_data)
        self.assertIsInstance(process_data['running'], bool)

        self.assertIn('status', process_data)
        self.assertIsInstance(process_data['status'], int)

        self.assertIn('returncode', process_data)
        self.assertIsInstance(process_data['returncode'], int)

    def test_endpoint_process_results(self):
        response = requests.get(
            self.urlBase + '/processes?count=-1',
            timeout = 5
        )
        self.assertEqual(response.status_code,200)
        process_list = response.json()
        if not len(process_list['result']):
            self.start_basic_run()
            response = requests.get(
                self.urlBase + '/processes?count=-1',
                timeout = 5
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
        process_list['result'].sort(key = lambda x:x['id'])
        while process_list['result'][-1]['running']:
            time.sleep(1)
            response = requests.get(
                self.urlBase + '/processes?count=-1',
                timeout = 5
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
            process_list['result'].sort(key = lambda x:x['id'])
        response = requests.get(
            'http://localhost:8080'+process_list['result'][-1]['results_url'] + '?count=-1',
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        results = response.json()
        self.assertIsInstance(results, dict)
        self.assertIsInstance(results['result'], list)
        for item in results['result']:
            self.assertIsInstance(item, dict)

            self.assertIn('description', item)
            self.assertIsInstance(item['description'], str)
            self.assertTrue(item['description'])

            self.assertIn('display_name', item)
            self.assertIsInstance(item['display_name'], str)
            self.assertTrue(item['display_name'])

            self.assertIn('fileID', item)
            self.assertIsInstance(item['fileID'], int)
            self.assertGreaterEqual(int(item['fileID']), 0)

            self.assertIn('rows', item)
            self.assertIsInstance(item['rows'], int)
            self.assertGreaterEqual(item['rows'], -1)

            self.assertIn('size', item)
            self.assertIsInstance(item['size'], int)

            self.assertIn('url', item)
            self.assertIsInstance(item['url'], str)
            self.assertTrue(item['url'])

        for process in process_list['result']:
            response = requests.get(
                # %2B = '+', %2C = ','
                'http://localhost:8080' + process['results_url'] + '?type=final',
                timeout = 5
            )
            self.assertEqual(response.status_code,200)
            results = response.json()
            for item in results['result']:
                self.assertTrue(re.search('final.tsv$', item['display_name']))

    def test_endpoint_process_results_data(self):
        response = requests.get(
            self.urlBase + '/processes',
            timeout = 5
        )
        self.assertEqual(response.status_code,200)
        process_list = response.json()
        if not len(process_list['result']):
            self.start_basic_run()
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
        process_list['result'].sort(key = lambda x:x['id'])
        while process_list['result'][-1]['running']:
            time.sleep(1)
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
            process_list['result'].sort(key = lambda x:x['id'])
        response = requests.get(
            'http://localhost:8080'+process_list['result'][-1]['results_url'],
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        results = response.json()
        for item in results['result']:
            if item['display_name'].endswith('.tsv') and item['rows']>0:
                response = requests.get(
                    'http://localhost:8080'+item['url'],
                    timeout=5,
                )
                self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
                data = response.json()
                self.assertIsInstance(data, list)
                for row in data:
                    self.assertIsInstance(row, dict)
                    self.assertIn('rowid', row)
                #check the cols endpoint
                response = requests.get(
                    'http://localhost:8080'+item['url']+'/cols',
                    timeout = 5
                )
                self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
                #check the schema endpoint
                response = requests.get(
                    'http://localhost:8080'+item['url']+'/schema',
                    timeout = 5
                )
                self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())

    def test_endpoint_stop(self):
        processID = self.start_basic_run()['processid']
        response = requests.get(
            self.urlBase+'/stop/%d'%processID,
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        response = requests.get(
            self.urlBase+'/processes/%d'%processID
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        data = response.json()
        self.assertIsInstance(data, dict)
        self.assertIn('running', data)
        self.assertFalse(data['running'])

    def test_endpoint_archive(self):
        processID = self.start_basic_run()['processid']
        time.sleep(1)
        response = requests.get(
            self.urlBase+'/processes/%d'%processID,
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        process_data = response.json()
        self.assertIsInstance(process_data, dict)
        self.assertIn('running', process_data)
        while process_data['running']:
            time.sleep(5)
            response = requests.get(
                self.urlBase+'/processes/%d'%processID,
                timeout = 5
            )
            self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
            process_data = response.json()
        response = requests.get(
            self.urlBase+'/archive/%d'%processID,
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        self.assertIsInstance(response.json(), str)

    def test_full_api_pipeline(self):
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
                'tdna_vaf':'40',
                'binding_threshold':3000,
                'force':True

            }
        )
        self.assertEqual(response.status_code, 201, response.url+' : '+response.content.decode())
        processID = response.json()['processid']
        self.assertIsInstance(processID, int)
        response = requests.get(
            self.urlBase+'/processes/%d'%processID,
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        process_data = response.json()
        self.assertIsInstance(process_data, dict)
        self.assertIn('running', process_data)
        while process_data['running']:
            time.sleep(5)
            response = requests.get(
                self.urlBase+'/processes/%d'%processID,
                timeout = 5
            )
            self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
            process_data = response.json()
        # time.sleep(1)
        # response = requests.get(
        #     self.urlBase+'/processes/%d'%processID,
        #     timeout = 5
        # )
        # self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        # process_data = response.json()
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
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        content = response.json()
        self.assertIsInstance(content, list)
        self.assertTrue(content)
        response = requests.get(
            'http://localhost:8080'+finaltsv['url']+'/cols',
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        mapping = response.json()
        self.assertIsInstance(mapping, dict)
        mapping['rowid'] = 'rowid'
        testlines = [row for row in reader]
        self.assertEqual(len(testlines), len(content), "Line count mismatch")
        for (testrow, outputrow) in zip(testlines, content):
            self.assertEqual({mapping[key] for key in outputrow}-testrow.keys(), {'rowid'})
            del outputrow['rowid']
            for key in outputrow:
                self.assertEqual(outputrow[key], parsedata(testrow[mapping[key]]), "Mismatch: %s"%key)
        raw_reader.close()

    def test_endpoint_allele(self):
        response = requests.get(
            self.urlBase+'/checkallele',
            timeout=5,
            params={
                'allele':'H2-IAb'
            }
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        self.assertTrue(response.json())
        response = requests.get(
            self.urlBase+'/checkallele',
            timeout=5,
            params={
                'allele':'NotARealAllele'
            }
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        self.assertFalse(response.json())

    def test_endpoint_validallele(self):
        response = requests.get(
            self.urlBase+'/validalleles',
            timeout=5,
            params={
                'prediction_algorithms':'NetMHC'
            }
        )
        self.assertEqual(response.status_code, 200)
        results = response.json()
        self.assertIsInstance(results, dict)
        self.assertTrue('NetMHC' in results)
        self.assertTrue(len(results['NetMHC']))

    def test_endpoint_validalgorithms(self):
        response = requests.get(
            self.urlBase+'/validalgorithms',
            timeout=5
        )
        self.assertEqual(response.status_code, 200)
        algorithms = response.json()
        self.assertIsInstance(algorithms, list)
        self.assertTrue('NetMHC' in algorithms)

    def test_duplicate_check(self):
        processID = self.start_basic_run()['processid']
        time.sleep(1)
        response = requests.get(
            self.urlBase+'/processes/%d'%processID,
            timeout = 5
        )
        self.assertEqual(response.status_code, 200)
        while response.json()['running']:
            time.sleep(2)
            response = requests.get(
                self.urlBase+'/processes/%d'%processID,
                timeout = 5
            )
            self.assertEqual(response.status_code, 200)
        response = requests.post(
            self.urlBase+'/staging',
            timeout = 5,
            data={
                'input':os.path.join(
                    self.test_data_directory,
                    'input.vcf'
                ),
                'samplename':'basic_run',
                'alleles':'HLA-E*01:01',
                'prediction_algorithms':'NetMHC',
            }
        )
        self.assertEqual(response.status_code, 400)
        self.assertRegex(
            response.json()['message'],
            re.compile(r"The given parameters match process \d+")
        )

    def test_endpoint_restart(self):
        processID = self.start_basic_run()['processid']
        time.sleep(1)
        response = requests.get(
            self.urlBase+'/processes/%d'%processID,
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        old_data = response.json()
        self.assertIsInstance(old_data, dict)
        self.assertIn('running', old_data)
        while old_data['running']:
            time.sleep(5)
            response = requests.get(
                self.urlBase+'/processes/%d'%processID,
                timeout = 5
            )
            self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
            old_data = response.json()
        response = requests.get(
            self.urlBase+'/restart/%d'%processID,
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
        process_data = response.json()
        self.assertIsInstance(process_data, dict)

        self.assertIn('attached', process_data)
        self.assertTrue(process_data['attached'])

        self.assertIn('command', process_data)
        self.assertEqual(process_data['command'], old_data['command'])

        self.assertIn('id', process_data)
        self.assertEqual(process_data['id'], old_data['id'])

        self.assertIn('output', process_data)
        self.assertEqual(process_data['output'], old_data['output'])

        self.assertIn('parameters', process_data)
        self.assertDictEqual(process_data['parameters'], old_data['parameters'])

        self.assertIn('pid', process_data)
        self.assertNotEqual(process_data['pid'], old_data['pid'])

        self.assertIn('results_url', process_data)
        self.assertEqual(process_data['results_url'], old_data['results_url'])

    def test_endpoint_dropbox(self):
        shutil.copyfile(
            os.path.join(
                self.test_data_directory,
                'Test.final.tsv'
            ),
            os.path.expanduser(os.path.join(
                '~',
                'pVAC-Seq',
                'dropbox',
                'Test.final.tsv'
            ))
        )
        time.sleep(1)
        response = requests.get(
            self.urlBase+'/dropbox',
            timeout = 5
        )
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())
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
            self.assertIsInstance(item['size'], int)

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
        self.assertEqual(response.status_code, 200, response.url+' : '+response.content.decode())

    # Some simple testing of filtering, sorting, and paging.
    # In the future, maybe add more thorough testing of these
    def test_sort(self):
        response = requests.get(
            self.urlBase + '/processes',
            timeout = 5,
            params={
                'count':-1
            }
        )
        self.assertEqual(response.status_code,200)
        process_list = response.json()
        if not len(process_list['result']):
            self.start_basic_run()
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5,
                params={
                    'count':-1
                }
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
        while process_list['result'][-1]['running']:
            time.sleep(1)
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5,
                params={
                    'count':-1
                }
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()

        for process in process_list['result']:
            response = requests.get(
                'http://localhost:8080' + process['results_url'],
                timeout = 5,
                params={
                    'sorting':'+size,+fileID',
                    'count':-1
                }
            )
            self.assertEqual(response.status_code,200)
            result_list = response.json()

            response2 = requests.get(
                'http://localhost:8080' + process['results_url'],
                timeout = 5,
                params={
                    'count':-1
                }
            )
            self.assertEqual(response2.status_code,200)
            normal_result_list = response2.json()
            normal_result_list['result'].sort(key = lambda x: (x['size'], x['fileID']))
            self.assertEqual(result_list['result'], normal_result_list['result'])

    def test_filter(self):
        response = requests.get(
            self.urlBase + '/processes',
            timeout = 5,
            params={
                'count':-1
            }
        )
        self.assertEqual(response.status_code,200)
        process_list = response.json()
        if not len(process_list['result']):
            self.start_basic_run()
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5,
                params={
                    'count':-1
                }
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
        while process_list['result'][-1]['running']:
            time.sleep(1)
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5,
                params={
                    'count':-1
                }
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()

        for process in process_list['result']:
            response = requests.get(
                'http://localhost:8080' + process['results_url'],
                timeout = 5,
                params={
                    'filters':'fileID>2'#,
                    #'count':-1
                }
            )
            filtered_results = response.json()
            self.assertEqual(response.status_code,200)
            for item in filtered_results['result']:
                self.assertGreater(item['fileID'], 2)

    def test_pagination(self):
        response = requests.get(
            self.urlBase + '/processes',
            timeout = 5,
            params={
                'count':-1
            }
        )
        self.assertEqual(response.status_code,200)
        process_list = response.json()
        if not len(process_list['result']):
            self.start_basic_run()
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5,
                params={
                    'count':-1
                }
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()
        while process_list['result'][-1]['running']:
            time.sleep(1)
            response = requests.get(
                self.urlBase + '/processes',
                timeout = 5,
                params={
                    'count':-1
                }
            )
            self.assertEqual(response.status_code,200)
            process_list = response.json()

        for process in process_list['result']:
            response = requests.get(
                'http://localhost:8080' + process['results_url'],
                timeout = 5,
                params={
                    'count':-1
                }
            )
            self.assertEqual(response.status_code,200)
            full_result_list = response.json()
            response = requests.get(
                'http://localhost:8080' + process['results_url'],
                timeout = 5,
                params={
                    'count':3,
                    'page':2
                }
            )
            self.assertEqual(response.status_code,200)
            paged_result = response.json()

            self.assertEqual(len(paged_result['result']), 3 if len(full_result_list['result']) >= 6 else len(full_result_list['result']) - 3)
            self.assertEqual(paged_result['result'][0] if len(paged_result['result']) else [], full_result_list['result'][3] if len(full_result_list['result']) > 3 else [])
            self.assertEqual(paged_result['_meta']['per_page'], 3)
            self.assertEqual(paged_result['_meta']['current_page'], 2)
            self.assertEqual(paged_result['_meta']['total_count'], len(full_result_list['result']))


