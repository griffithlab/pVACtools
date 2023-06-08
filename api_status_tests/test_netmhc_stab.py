import unittest
import os
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pvactools.lib.netmhc_stab import NetMHCStab
from tests.utils import *


class NetMHCstabpanApiStatusTest(unittest.TestCase):
    def test_netmhc_stab_api_runs_and_returns_expected_value(self):
        staging_file = tempfile.NamedTemporaryFile(mode='w+')
        record = SeqRecord(Seq("ASTPGHTIIYEAVCLHNDRTTIP"), id="1", description="1")
        SeqIO.write([record], staging_file.name, "fasta")
        staging_file.seek(0)
        response = NetMHCStab.query_netmhcstabpan_server(None, staging_file, 9, "HLA-A02:01")
        expected_file_name = os.path.join(pvactools_directory(), 'api_status_tests', 'netmhcstab_response.html')
        with open(expected_file_name, 'r') as expected_fh:
            actual_content = response.content.decode('ascii')
            expected_content = expected_fh.read()
            self.assertEqual(expected_content, actual_content)
