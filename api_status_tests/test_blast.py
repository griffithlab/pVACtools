import unittest
from Bio.Blast import NCBIWWW

from tests.utils import *


class BlastApiStatusTest(unittest.TestCase):
    def test_blast_api_runs_and_returns_expected_value(self):
        result_handle = NCBIWWW.qblast("blastp", "refseq_select_prot", "AAIMYVPALGWEFLAFTRLTSELNFLLQEID", entrez_query="Homo Sapiens [Organism]", word_size=7, gapcosts='32767 32767')
        expected_file_name = os.path.join(pvactools_directory(), 'api_status_tests', 'blast_response.xml')
        with open(expected_file_name, 'r') as expected_fh:
            actual_content = result_handle.read()
            expected_content = expected_fh.read()
            actual_content_without_query_id = '\n'.join(line for line in actual_content.split('\n') if 'BlastOutput_query-ID' not in line and 'Iteration_query-ID' not in line)
            expected_content_without_query_id = '\n'.join(line for line in expected_content.split('\n') if 'BlastOutput_query-ID' not in line and 'Iteration_query-ID' not in line)
            self.assertEqual(expected_content_without_query_id, actual_content_without_query_id)
