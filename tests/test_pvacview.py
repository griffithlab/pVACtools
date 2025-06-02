import unittest
import unittest.mock
import os
import re
import sys
import time
import tempfile
import py_compile
import subprocess
from subprocess import PIPE
from subprocess import run as subprocess_run
import socket
import argparse

from pvactools.tools.pvacview import *
import pvactools.tools.pvacview.main as pvacview_main
from tests.utils import *

class PvacviewTests(unittest.TestCase):
    def test_pvacview_compiles(self):
        compiled_pvac_path = py_compile.compile(os.path.join(
            pvactools_directory(),
            'pvactools',
            'tools',
            'pvacview',
            "main.py"
        ))
        self.assertTrue(compiled_pvac_path)

    def test_parser(self):
        parser = pvacview_main.define_parser()
        self.assertEqual(type(parser), argparse.ArgumentParser)

    def test_pvacview_commands(self):
        pvac_script_path = os.path.join(
            pvactools_directory(),
            'pvactools',
            'tools',
            'pvacview',
            "main.py"
            )
        usage_search = re.compile(r"usage: ")
        result = subprocess_run([
            sys.executable,
            pvac_script_path,
            "run",
            '-h'
        ], shell=False, stdout=PIPE)
        self.assertFalse(result.returncode, "Failed `pvacseq run -h`")
        self.assertRegex(result.stdout.decode(), usage_search)

    def test_run_compiles(self):
        compiled_run_path = py_compile.compile(os.path.join(
            pvactools_directory(),
            'pvactools',
            "tools",
            "pvacview",
            "run.py"
        ))
        self.assertTrue(compiled_run_path)

    @unittest.skipIf(os.getenv('GITHUB_ACTIONS') == 'true', "Skipping test on GitHub Actions")
    def test_run(self):
        pvac_script_path = os.path.join(
            pvactools_directory(),
            'pvactools',
            'tools',
            'pvacview',
            "main.py"
            )
        example_data_path = os.path.join(pvactools_directory(), "pvactools", "tools", "pvacview")
        server_process = subprocess.Popen(
            [sys.executable, pvac_script_path, "run", example_data_path],
            stdout=PIPE,
            stderr=PIPE
        )
        time.sleep(10)
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.connect(('localhost', 3333))
            s.sendall(b"GET / HTTP/1.1\r\nHost: localhost\r\n\r\n")
            data = s.recv(1024)
            s.close()
            self.assertIn(b"200 OK", data, "Server is not responding correctly")
        except socket.error:
            self.fail("Server did not start.")
            self.fail(f"Error communicating with server: {e}")
        server_process.terminate()
        server_process.wait()
