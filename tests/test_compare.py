import unittest
import subprocess
import requests
import time
import py_compile
import argparse
import os
import pvactools.tools.compare as pvaccompare_main

# python -m unittest tests/test_compare.py
# python -m unittest discover -s tests
class TestRunCompare(unittest.TestCase):
    def setUp(self):
        self.server_process = subprocess.Popen(
            ["python", "-m", "pvactools.tools.pvaccompare.server"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        time.sleep(1)

    def tearDown(self):
        self.server_process.terminate()
        self.server_process.wait()

    def test_compare_server(self):
        url = "http://localhost:8080/pvactools/tools/pvaccompare/html_report/main.html"
        response = requests.get(url)
        self.assertEqual(response.status_code, 200)

    def test_compare_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            'pvactools',
            'tools',
            "compare.py"
        ))
        self.assertTrue(compiled_pvac_path)
    
    def test_compare_parser(self):
        parser = pvaccompare_main.define_parser()
        self.assertEqual(type(parser), argparse.ArgumentParser)
