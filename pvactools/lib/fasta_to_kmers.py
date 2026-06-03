import os
import re
import pandas as pd
import numpy as np
from pyfaidx import Fasta
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, translate

from pvactools.lib.run_utils import *

class FastaToKmers:
    def __init__(self, **kwargs):
        self.tscript_fasta = Fasta(kwargs['fasta'])
        self.fasta_path    = kwargs['fasta']
        self.output_dir    = kwargs['output_dir']
        self.epitope_length = kwargs['epitope_length']
        self.sample_name   = kwargs['sample_name']

    def create_kmers(self, seq_name):
        supported_aas = supported_amino_acids()
        kmer_dict = {}
        # using personalized fasta
        sequence = str(self.tscript_fasta[seq_name])
        # loop over entire sequence (doing this for WT and ALT) i == position in peptide
        for i in range(len(sequence)):
            # grab kmer sequence
            k = sequence[i:self.epitope_length+i]
            if len(k) < self.epitope_length:
                continue
            #kmer contains unsupported amino acids
            if not all([c in supported_aas for c in k]):
                logging.warning("Record {} contains unsupported amino acids. Skipping.".format(k))
                continue
            # add entry to dictionary
            kmer_dict[i] = k

        return kmer_dict

    def create_kmer_dict(self, splice_site_name):
        # wt and mut each get kmer dict
        wt_dict = self.create_kmers(f'WT.{splice_site_name}')
        mut_dict = self.create_kmers(f'ALT.{splice_site_name}')
        splice_type = splice_site_name.rsplit('.', 1)[1]
        final_kmers = {}
        for i in range(len(mut_dict)):
            mt_seq = mut_dict[i]
            if i < len(wt_dict):
                wt_seq = wt_dict[i]
            else:
                wt_seq = ""
            min_match = min_match_count(len(mt_seq))
            diff = len(wt_dict) - len(mut_dict)
            alt_i = i + diff
            if alt_i < 0:
                alt_wt_seq = ""
            else:
                alt_wt_seq = wt_dict[i+diff]
            #Skip MT kmer if it occurs in the WT dict
            if mt_seq in wt_dict.values():
                continue
            else:
                left_match_count = determine_consecutive_matches_from_left(mt_seq, wt_seq)
                if splice_type == 'inframe_splice_site':
                    right_match_count = determine_consecutive_matches_from_right(mt_seq, alt_wt_seq)
                    if left_match_count >= right_match_count:
                        wt_seq_to_consider = wt_seq
                    else:
                        wt_seq_to_consider = alt_wt_seq
                    total_match_count = determine_total_matches(mt_seq, wt_seq_to_consider)
                    if total_match_count >= min_match:
                        final_wt_seq = wt_seq_to_consider
                    else:
                        final_wt_seq = None
                elif splice_type == 'frameshift_splice_site':
                    total_match_count = determine_total_matches(mt_seq, wt_seq)
                    if total_match_count >= min_match:
                        final_wt_seq = wt_seq
                    else:
                        final_wt_seq = None
                final_kmers[f'ALT.{splice_site_name}|{i}'] = mt_seq
                if final_wt_seq is not None:
                    final_kmers[f'WT.{splice_site_name}|{i}'] = final_wt_seq
        return final_kmers

    def prefix(self):
        return "ALT."

    def loop_through_tscripts(self):
        all_kmers = {}
        # all MT fasta headers
        alt_fasta_keys = [k for k in self.tscript_fasta.keys() if k.startswith(self.prefix())]
        for alt_name in alt_fasta_keys:
            splice_site_name = alt_name.removeprefix(self.prefix())
            # get final mutated kmer list from save_kmer_dicts()
            all_kmers.update(self.create_kmer_dict(splice_site_name))
        return all_kmers


    def create_epitope_fasta(self, all_kmers):
        records = []
        for index, seq in all_kmers.items():
            records.append(SeqRecord(Seq(seq), id=index, description=""))
        output_file = f'{self.output_dir}/{self.sample_name}.{self.epitope_length}.fa'
        SeqIO.write(records, output_file, "fasta")

    def execute(self):
        all_kmers = self.loop_through_tscripts()
        if len(all_kmers) > 0:
            self.create_epitope_fasta(all_kmers)

class FusionFastaToKmers(FastaToKmers):
    def prefix(self):
        return "MT."

    def create_kmer_dict(self, index):
        # wt and mut each get kmer dict
        five_wt_dict = self.create_kmers(f'WT5.{index}')
        fusion_type = index.rsplit('.', 2)[1]
        if fusion_type == 'inframe_fusion':
            three_wt_dict = self.create_kmers(f'WT3.{index}')
        mut_dict = self.create_kmers(f'MT.{index}')

        final_kmers = {}
        for i in range(len(mut_dict)):
            mt_seq = mut_dict[i]
            if mt_seq in five_wt_dict.values():
                continue
            if fusion_type == 'inframe_fusion' and mt_seq in three_wt_dict.values():
                continue

            min_match = min_match_count(len(mt_seq))
            if fusion_type == 'frameshift_fusion':
                if i in five_wt_dict:
                    five_wt_seq = five_wt_dict[i]
                    total_match_count = determine_total_matches(mt_seq, five_wt_seq)
                    if total_match_count >= min_match:
                        final_wt_seq = five_wt_seq
                    else:
                        final_wt_seq = None
                else:
                    final_wt_seq = None
            else:
                if i in five_wt_dict:
                    five_wt_seq = five_wt_dict[i]
                    left_match_count = determine_consecutive_matches_from_left(mt_seq, five_wt_seq)
                else:
                    left_match_count = 0
                alt_i = list(three_wt_dict.keys())[-1] - list(mut_dict.keys())[-1] + i
                if alt_i in three_wt_dict:
                    three_wt_seq = three_wt_dict[alt_i]
                    right_match_count = determine_consecutive_matches_from_right(mt_seq, three_wt_seq)
                else:
                    right_match_count = 0
                if left_match_count >= right_match_count:
                    wt_seq_to_consider = five_wt_seq
                else:
                    wt_seq_to_consider = three_wt_seq
                total_match_count = determine_total_matches(mt_seq, wt_seq_to_consider)
                if total_match_count >= min_match:
                    final_wt_seq = wt_seq_to_consider
                else:
                    final_wt_seq = None

            final_kmers[f'MT.{index}|{i}'] = mt_seq
            if final_wt_seq is not None:
                final_kmers[f'WT.{index}|{i}'] = final_wt_seq

        return final_kmers

