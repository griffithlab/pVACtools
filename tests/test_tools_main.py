import unittest
import tools.main as tools_main
import tools.download_cwls as download_cwls
import argparse
import tempfile
import os

class ToolsMainTests(unittest.TestCase):
    def test_parser(self):
        parser = tools_main.define_parser()
        self.assertEqual(type(parser), argparse.ArgumentParser)

    def test_download_cwls(self):
        output_dir = tempfile.TemporaryDirectory()
        download_cwls.main([output_dir.name])
        self.assertTrue(os.path.exists(os.path.join(output_dir.name, 'pvactools_cwls', 'pvacseq.cwl')))
        self.assertTrue(os.path.exists(os.path.join(output_dir.name, 'pvactools_cwls', 'pvacfuse.cwl')))
        self.assertTrue(os.path.exists(os.path.join(output_dir.name, 'pvactools_cwls', 'pvacvector.cwl')))
