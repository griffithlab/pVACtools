import argparse
import sys
import requests
from requests.exceptions import Timeout
import csv
import tempfile
import re
import os
from time import sleep
import random
import logging
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import pvactools.lib.run_utils
from pvactools.lib.prediction_class import MHCI
import pvactools.lib.sort

class NetMHCStab:
    def __init__(self, input_file, output_file, file_type='pVACseq', top_score_metric='median'):
        self.input_file = input_file
        self.output_file = output_file
        if file_type == 'pVACseq':
            self.epitope_seq_column_name = 'MT Epitope Seq'
        else:
            self.epitope_seq_column_name = 'Epitope Seq'
        self.file_type = file_type
        self.top_score_metric = top_score_metric

    def execute(self):
        mhci_alleles = MHCI.all_valid_allele_names()
        observed_alleles = self.observed_alleles()
        alleles = self.valid_alleles(list(set(observed_alleles).intersection(set(mhci_alleles))))
        invalid_alleles = list(set(observed_alleles) - set(alleles))
        lengths = self.observed_epitope_lengths()

        result_delimiter = re.compile(r'-{20,}')
        fail_searcher = re.compile(r'(Failed run|Problematic input:|Configuration error)')
        rejected_searcher = re.compile(r'status: rejected')
        success_searcher = re.compile(r'Rank Threshold for Strong binding peptides')
        cannot_open_file_searcher = re.compile(r'Cannot open file')
        allele_searcher = re.compile(r'^(.*?) : Distance to trai?ning data\s+(\d.\d+).*? nearest neighbor (.*?)\)$', re.MULTILINE)
        with open(self.output_file, 'w') as output_fh:
            headers = pd.read_csv(self.input_file, delimiter="\t", nrows=0).columns.tolist()
            writer = csv.DictWriter(
                output_fh,
                headers+['Predicted Stability', 'Half Life', 'Stability Rank', 'NetMHCstab allele'],
                delimiter='\t',
                lineterminator='\n'
            )
            writer.writeheader()

            output_lines = []
            for allele in alleles:
                for length in lengths:
                    df = (pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False, na_values="NA", keep_default_na=False)
                            [lambda x: (x['HLA Allele'] == allele) & (x[self.epitope_seq_column_name].str.len() == length) ])
                    df = df.fillna('NA')
                    if len(df) == 0:
                        continue

                    netmhcstabpan_allele = allele.replace('*', '')
                    x = 0
                    chunk_count = int(len(df)/100)
                    if chunk_count == 0:
                        chunk_count = 1
                    for chunk in np.array_split(df, chunk_count):
                        staging_file = tempfile.NamedTemporaryFile(mode='w+')
                        current_buffer = {}
                        records = []
                        for index, line in chunk.iterrows():
                            sequence_id = ('%010x'%x)[-10:]
                            record_new = SeqRecord(Seq(line[self.epitope_seq_column_name], IUPAC.protein), id=sequence_id, description=sequence_id)
                            records.append(record_new)
                            current_buffer[sequence_id] = line.to_dict()
                            x+=1
                        SeqIO.write(records, staging_file.name, "fasta")
                        staging_file.seek(0)
                        response = self.query_netmhcstabpan_server(staging_file, length, netmhcstabpan_allele)

                        if fail_searcher.search(response.content.decode()):
                            raise Exception("NetMHCstabpan encountered an error during processing.\n{}".format(response.content.decode()))

                        while rejected_searcher.search(response.content.decode()):
                            logging.warning("Too many jobs submitted to NetMHCstabpan server. Waiting to retry.")
                            sleep(random.randint(5, 10))
                            staging_file.seek(0)
                            response = self.query_netmhcstabpan_server(staging_file, length, netmhcstabpan_allele)

                        if cannot_open_file_searcher.search(response.content.decode()):
                            sleep(random.randint(5, 10))
                            staging_file.seek(0)
                            response = self.query_netmhcstabpan_server(staging_file, length, netmhcstabpan_allele)
                            while rejected_searcher.search(response.content.decode()):
                                sleep(random.randint(5, 10))
                                staging_file.seek(0)
                                response = self.query_netmhcstabpan_server(staging_file, length, netmhcstabpan_allele)
                            if cannot_open_file_searcher.search(response.content.decode()):
                                raise Exception("NetMHCstabpan server was unable to read the submitted fasta file:\n{}.".format(staging_file.read()))

                        if success_searcher.search(response.content.decode()):
                            allele_map = {}
                            for item in allele_searcher.findall(response.content.decode()):
                                allele_map[item[0]] = "{} (distance: {})".format(item[2], item[1])
                                if item[1] != "0.000":
                                    print("NetMHCstabpan substituted {} for {} (distance: {})".format(item[2], item[0], item[1]))
                            results = [item.strip() for item in result_delimiter.split(response.content.decode())]
                            if len(results) == 0:
                                raise Exception("Unexpected return value from NetMHCstabpan server. Unable to parse response.\n{}".format(response.content.decode()))
                            data_for_sequence_id = {}
                            for i in range(2, len(results), 4): #examine only the parts we want, skipping all else
                                for result_line in results[i].split('\n'):
                                    data = [word for word in result_line.strip().split(' ') if len(word)]
                                    data_for_sequence_id[data[3]] = data

                            for sequence_id, line in current_buffer.items():
                                data = data_for_sequence_id[sequence_id]
                                line.update({
                                    'Predicted Stability':data[4],
                                    'Half Life':data[5],
                                    'Stability Rank':data[6],
                                    'NetMHCstab allele':allele_map[netmhcstabpan_allele]
                                })
                                output_lines.append(line)

                        else:
                            raise Exception("Unexpected return value from NetMHCstabpan server. Unable to parse response.\n{}".format(response.content.decode()))

            for allele in invalid_alleles:
                df = (pd.read_csv(self.input_file, delimiter='\t', float_precision='high', low_memory=False, na_values="NA", keep_default_na=False)
                        [lambda x: (x['HLA Allele'] == allele) ])
                df['Predicted Stability'] = 'NA'
                df['Half Life'] = 'NA'
                df['Stability Rank'] = 'NA'
                df['NetMHCstab allele'] = 'NA'
                output_lines.extend(df.to_dict('records'))

            if self.file_type == 'pVACseq':
                sorted_lines = pvactools.lib.sort.default_sort_from_pd_dict(output_lines, self.top_score_metric)
            else:
                sorted_lines = pvactools.lib.sort.pvacbind_sort(output_lines, self.top_score_metric)
            writer.writerows(sorted_lines)

    def query_netmhcstabpan_server(self, staging_file, peptide_length, allele):
        try:
            response = requests.post(
                "https://services.healthtech.dtu.dk/cgi-bin/webface2.cgi",
                files={'SEQSUB':(staging_file.name, staging_file, 'text/plain')},
                data = {
                    'configfile':'/var/www/html/services/NetMHCstabpan-1.0/webface.cf',
                    'inp':'0',
                    'len': peptide_length,
                    'master':'1',
                    'slave0':allele,
                    'allele':allele,
                    'thrs':'0.5',
                    'thrw': '2',
                    'incaff': '0',
                    'sort1':'-1',
                    'waff':'0.8',
                    'sort2':'-1',
                },
                timeout=(10,60)
            )
        except Timeout:
            raise Exception("Timeout while posting request to NetMHCstabpan server. The server may be unresponsive. Please try again later.")
        if response.status_code != 200:
            raise Exception("Error posting request to NetMHCstabpan server.\n{}".format(response.content.decode()))

        jobid_searcher = re.compile(r'<!-- jobid: [0-9a-fA-F]*? status: (queued|active)')
        while jobid_searcher.search(response.content.decode()):
            sleep(10)
            try:
                response = requests.get(response.url, timeout=(10,60))
            except Timeout:
                raise Exception("Timeout while posting request to NetMHCstabpan server. The server may be unresponsive. Please try again later.")
            if response.status_code != 200:
                raise Exception("Error posting request to NetMHCstabpan server.\n{}".format(response.content.decode()))
        return response

    def valid_alleles(self, alleles):
        invalid_searcher = re.compile(r'cannot be found in hla_pseudo list')
        valid_alleles = []
        for allele in alleles:
            staging_file = tempfile.NamedTemporaryFile(mode='w+')
            records = [SeqRecord(Seq("ASTPGHTIIYEAVCLHNDRTTIP", IUPAC.protein), id="0", description="0")]
            SeqIO.write(records, staging_file.name, "fasta")
            staging_file.seek(0)
            response = self.query_netmhcstabpan_server(staging_file, 9, allele.replace("*", ""))

            if not invalid_searcher.search(response.content.decode()):
                valid_alleles.append(allele)
        return valid_alleles

    def observed_alleles(self):
        return np.sort(pd.read_csv(self.input_file, delimiter="\t", usecols=["HLA Allele"])['HLA Allele'].unique())[::-1]

    def observed_epitope_lengths(self):
        epitopes = pd.read_csv(self.input_file, delimiter="\t", usecols=[self.epitope_seq_column_name])[self.epitope_seq_column_name].unique()
        return list(set([len(e) for e in epitopes]))

    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            "%s netmhc_stab" % tool,
            description="Add stability predictions to predicted neoepitopes.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="Input filtered file with predicted epitopes."
        )
        parser.add_argument(
            'output_file',
            help="Output TSV filename for putative neoepitopes."
        )
        parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            default='median',
            help="The ic50 scoring metric to use when sorting epitopes. "
                 + "lowest: Use the best MT Score and Corresponding Fold Change (i.e. the lowest MT ic50 binding score and corresponding fold change of all chosen prediction methods). "
                 + "median: Use the median MT Score and Median Fold Change (i.e. the  median MT ic50 binding score and fold change of all chosen prediction methods)."
        )
        return parser
