import unittest
import os
import sys
import tempfile
from filecmp import cmp
import py_compile
import pandas as pd
from pyfaidx import FastaVariant

from pvactools.lib.junction_to_fasta import JunctionToFasta
from tests.utils import *

#python -m unittest tests/test_pvacsplice_junction_to_fasta.py
class JunctionToFastaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.python = sys.executable
        cls.executable = os.path.join(pvactools_directory(), "pvactools", "tools", "pvacsplice", "junction_to_fasta.py")
        cls.test_data_dir = os.path.join(pvactools_directory(), "tests", "test_data", "pvacsplice")

    def module_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    # using GBM data instead? - how about i run one with current gtf and fa - one of every type going through the pipeline
    def test_junction_to_fasta_runs_and_produces_expected_output(self):
        combined_df   = pd.read_csv(os.path.join(self.test_data_dir, 'results', 'Test.combined.tsv'), sep='\t', na_values="NA", dtype={'transcript_support_level': str})
        gtf_df        = pd.read_csv(os.path.join(self.test_data_dir, 'results', 'Test.gtf.tsv'), sep='\t', na_values="NA", dtype={'transcript_support_level': str})
        annotated_vcf = os.path.join(pvactools_directory(), self.test_data_dir, 'inputs', 'annotated.expression_chr1.vcf.gz')
        sample_name   = 'HCC1395_TUMOR_DNA'
        fasta_path = os.path.join(self.test_data_dir, 'inputs', 'all_sequences_chr1.fa.gz')
        unzipped_fasta_file = gunzip_file(fasta_path)
        output_dir  = tempfile.TemporaryDirectory() #output_dir.name
        output_file = os.path.join(output_dir.name, 'sample_transcripts.fa')

        for i in combined_df.index.to_list():
            print(i)
            junction = combined_df.loc[[i], :]
            for row in junction.itertuples():
                junction_params = {
                    'fasta_path'     : unzipped_fasta_file,
                    'junction_df'    : combined_df,
                    'gtf_df'         : gtf_df,
                    'tscript_id'     : row.transcript_id,
                    'chrom'          : row.junction_chrom,
                    'junction_name'  : row.name,
                    'junction_coors' : [row.junction_start, row.junction_stop],
                    'fasta_index'    : row.index,
                    'variant_info'   : row.variant_info,
                    'anchor'         : row.anchor,
                    'strand'         : row.strand,
                    'gene_name'      : row.gene_name,
                    'output_file'    : output_file,
                    'output_dir'     : output_dir.name,
                    'sample_name'    : sample_name,
                    'vcf'            : annotated_vcf,
                }
                junctions = JunctionToFasta(**junction_params)
                wt = junctions.create_wt_df()
                if wt.empty:
                    continue
                alt = junctions.create_alt_df()
                if alt.empty:
                    continue
                wt_aa, wt_fs = junctions.get_aa_sequence(wt)
                alt_aa, alt_fs = junctions.get_aa_sequence(alt)
                if wt_aa == '' or alt_aa == '':
                    print('No amino acid sequence was produced. Skipping.')
                    continue
                # creates output transcript fasta
                junctions.create_sequence_fasta(wt_aa, alt_aa)

        expected_file = os.path.join(self.test_data_dir, 'results', 'Test.transcripts.fa')

        self.assertTrue(cmp(
                output_file,
                expected_file),
                "files don't match {} - {}".format(output_file, expected_file)
            )

        os.unlink(unzipped_fasta_file)

        output_dir.cleanup()
