import unittest
import os
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pvactools.lib.net_chop import NetChop
from tests.utils import *


class NetChopApiStatusTest(unittest.TestCase):
    def test_netchop_api_runs_and_returns_expected_value(self):
        staging_file = tempfile.NamedTemporaryFile(mode='w+')
        record = SeqRecord(Seq("ASTPGHTIIYEAVCLHNDRTTIP"), id="1", description="1")
        SeqIO.write([record], staging_file.name, "fasta")
        staging_file.seek(0)
        http = NetChop.setup_adapter(None)
        response = NetChop(None, None, None).query_netchop_server(http, staging_file, '0', 0.5)
        expected_file_name = os.path.join(pvactools_directory(), 'api_status_tests', 'netchop_response.html')
        with open(expected_file_name, 'r') as expected_fh:
            actual_content = response.content.decode('ascii')
            expected_content = expected_fh.read()
            self.assertEqual(expected_content, actual_content)
