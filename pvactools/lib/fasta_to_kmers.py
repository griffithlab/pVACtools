import os
import re
import pandas as pd
import numpy as np
from pyfaidx import Fasta
import logging

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
            final_seq_name = f'{seq_name};{i+1}'
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
            wt_seq = wt_dict[i]
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

    def loop_through_tscripts(self):
        unique_kmers = {}
        # all MT fasta headers
        alt_fasta_keys = [k for k in self.tscript_fasta.keys() if k.startswith('ALT.')]
        for alt_name in alt_fasta_keys:
            splice_site_name = alt_name.removeprefix("ALT.")
            # get final mutated kmer list from save_kmer_dicts()
            final_kmers = self.create_kmer_dict(splice_site_name)
            if not final_kmers:
                print(f'No unique kmers found for {splice_site_name}')
                continue
            # create master dict of unique kmers: index(es)
            for index, kmer in final_kmers.items():
                if kmer not in unique_kmers.keys():
                    unique_kmers[kmer] = [index]
                else:
                    unique_kmers[kmer].append(index)
        return unique_kmers


    def create_fasta_df(self, kmers_dict):
        # joined indexes for each unique kmer - to format for df since lists are dif sizes
        fasta_info = {k:','.join(sorted(v)) for k,v in kmers_dict.items()}
        # add header cols
        headers_dict = {'peptide': fasta_info.keys(), 'name': fasta_info.values()} # maybe sort here before turned into df
        # convert to df
        fasta_df = pd.DataFrame.from_dict(headers_dict) # peptide, name made into a df
        return fasta_df

    def create_epitope_fastas(self, fasta_df):
        fasta_df['name'] = fasta_df['name'].str.replace(';', '.')
        # sort the names
        len_subset = fasta_df.sort_values(by=['name'])
        # 1 file per kmer length
        output_file = f'{self.output_dir}/{self.sample_name}.{self.epitope_length}.fa'
        # loop over rows in subset df
        for row in len_subset.itertuples():
            # fasta entry
            write_str = f'>{row.name}\n{row.peptide}\n'
            # don't duplicate entries
            if os.path.exists(output_file):
                with open(output_file, "r+") as f:
                    dup_content = re.search(row.peptide, f.read())
                    if not dup_content:
                        f.write(write_str)
            else:
                with open(output_file, "w") as e:
                    e.write(write_str)

    def execute(self):
        unique_kmers = self.loop_through_tscripts()
        if len(unique_kmers) > 0:
            # key: peptide value: list of ids
            fasta_df = self.create_fasta_df(unique_kmers)
            self.create_epitope_fastas(fasta_df)
