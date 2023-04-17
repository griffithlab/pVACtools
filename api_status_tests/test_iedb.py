import unittest
import os
import tempfile
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tests.utils import *


class IEDBApiStatusTest(unittest.TestCase):
    def test_iedb_class_i_api_runs_and_returns_expected_value(self):
        staging_file = tempfile.NamedTemporaryFile(mode='w+')
        record = SeqRecord(Seq("ASTPGHTIIYEAVCLHNDRTTIP"), id="1", description="1")
        SeqIO.write([record], staging_file.name, "fasta")
        staging_file.seek(0)

        url = 'http://tools-cluster-interface.iedb.org/tools_api/mhci/'
        data = {
            'sequence_text': staging_file.read(),
            'method':        'ann',
            'allele':        'HLA-A*02:01',
            'length':        '9',
            'user_tool':     'pVac-seq',
        }

        response = requests.post(url, data=data)
        expected_file_name = os.path.join(pvactools_directory(), 'api_status_tests', 'iedb_class_i_response.tsv')
        with open(expected_file_name, 'r') as expected_fh:
            actual_content = response.content.decode('ascii')
            expected_content = expected_fh.read()
            self.assertEqual(expected_content, actual_content)

    def test_iedb_class_ii_api_runs_and_returns_expected_value(self):
        staging_file = tempfile.NamedTemporaryFile(mode='w+')
        record = SeqRecord(Seq("ASTPGHTIIYEAVCLHNDRTTIP"), id="1", description="1")
        SeqIO.write([record], staging_file.name, "fasta")
        staging_file.seek(0)

        url = 'http://tools-cluster-interface.iedb.org/tools_api/mhcii/'
        data = {
            'sequence_text': staging_file.read(),
            'method':        'nn_align',
            'allele':        'DPA1*01:03/DPB1*02:01',
            'length':        '15',
            'user_tool':     'pVac-seq',
        }

        response = requests.post(url, data=data)
        expected_file_name = os.path.join(pvactools_directory(), 'api_status_tests', 'iedb_class_ii_response.tsv')
        with open(expected_file_name, 'r') as expected_fh:
            actual_content = response.content.decode('ascii')
            expected_content = expected_fh.read()
            self.assertEqual(expected_content, actual_content)
