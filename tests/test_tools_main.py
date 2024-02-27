import unittest
import argparse
import tempfile
import os

import pvactools.tools.main as tools_main
import pvactools.tools.download_cwls as download_cwls
import pvactools.tools.download_wdls as download_wdls

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

    def test_download_wdls(self):
        output_dir = tempfile.TemporaryDirectory()
        download_wdls.main([output_dir.name])
        self.assertTrue(os.path.exists(os.path.join(output_dir.name, 'pvactools_wdls', 'pvacseq.wdl')))
        self.assertTrue(os.path.exists(os.path.join(output_dir.name, 'pvactools_wdls', 'pvacfuse.wdl')))
