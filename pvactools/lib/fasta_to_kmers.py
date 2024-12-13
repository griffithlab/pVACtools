import os
import re
import pandas as pd
import numpy as np
from pyfaidx import Fasta


class FastaToKmers:
    def __init__(self, **kwargs):
        self.tscript_fasta = Fasta(kwargs['fasta'])
        self.fasta_path    = kwargs['fasta']
        self.output_dir    = kwargs['output_dir']
        self.class_i_epitope_length  = kwargs['class_i_epitope_length']
        self.class_ii_epitope_length = kwargs['class_ii_epitope_length']
        self.class_i_hla   = kwargs['class_i_hla']
        self.class_ii_hla  = kwargs['class_ii_hla']
        self.sample_name   = kwargs['sample_name']
        self.final_lengths = self.choose_final_lengths()

    def create_kmers(self, seq_name):
        kmer_dict = {}
        # using personalized fasta
        sequence = str(self.tscript_fasta[seq_name])
        # loop over entire sequence (doing this for WT and ALT) i == position in peptide
        for i in range(len(sequence)):
            # loop over lengths
            for x in self.final_lengths:
                final_seq_name = f'{seq_name};{i+1}'
                # grab kmer sequence
                k = sequence[i:x+i]
                # if kmer matches target len
                # if k seq does not include X (any aa); continue
                match = re.search('X', k)
                if not match and len(k) == x:
                    # add entry to dictionary
                    kmer_dict[k] = final_seq_name

        return kmer_dict

    def create_unique_kmer_dict(self, wt_name, mt_name):
        # wt and mut each get kmer dict
        wt_dict = self.create_kmers(wt_name)
        mut_dict = self.create_kmers(mt_name)
        # all kmers not in wt kmers
        final_kmers = {k: v.replace('ALT.', '') for k, v in mut_dict.items() if k not in list(wt_dict.keys())}
        if len(final_kmers) == 0:
            final_kmers = {}
        return final_kmers

    def loop_through_tscripts(self):
        unique_kmers = {}
        # all fasta headers (WT and ALT)
        fasta_keys = list(self.tscript_fasta.keys())
        # take off WT. and ALT. prefixes from fasta_keys to de-duplicate
        unique_keys = sorted(set([x.split('.', 1)[1] for x in fasta_keys]))
        for key in unique_keys:
            # selecting WT and ALT for each seq pair (bc calling fasta_keys not unique_keys)
            wt_name, alt_name = [x for x in fasta_keys if key in x]
            # get final mutated kmer list from save_kmer_dicts()
            final_kmers = self.create_unique_kmer_dict(wt_name, alt_name)
            if not final_kmers:
                print(f'No unique kmers found for {key}')
                continue
            # create master dict of unique kmers: index(es)
            for k,v in final_kmers.items():
                if k not in unique_kmers.keys():
                    unique_kmers[k] = [v]
                else:
                    unique_kmers[k].append(v)
        return unique_kmers


    def create_index_file(self, kmers_dict):
        # joined indexes for each unique kmer - to format for df since lists are dif sizes
        fasta_info = {k:','.join(sorted(v)) for k,v in kmers_dict.items()}
        # add header cols
        headers_dict = {'peptide': fasta_info.keys(), 'name': fasta_info.values()} # maybe sort here before turned into df
        # convert to df
        fasta_df = pd.DataFrame.from_dict(headers_dict) # peptide, name made into a df
        # save to sep dictionary
        index_df = fasta_df.copy()
        # convert string to list and explode so each index is on a sep line
        index_df['name'] = index_df['name'].str.split(',')
        index_df = index_df.explode('name')
        # now split each index into index, pos, len
        index_df[['junction_index', 'transcript_position']] = index_df['name'].str.split(';', expand=True)
        index_df['name'] = index_df['name'].str.replace(';', '.')
        # create a tsv file

        return fasta_df

    def choose_final_lengths(self):
        if not self.class_i_hla:
            lengths = self.class_ii_epitope_length
        elif not self.class_ii_hla:
            lengths = self.class_i_epitope_length
        else:
            lengths = self.class_i_epitope_length + self.class_ii_epitope_length
        return lengths

    def create_epitope_fastas(self, fasta_df):
        fasta_df['name'] = fasta_df['name'].str.replace(';', '.')
        # sort the names
        # for only 1 length at a time
        for x in self.final_lengths:
            len_subset = fasta_df[fasta_df['peptide'].str.len() == x].sort_values(by=['name'])
            # 1 file per kmer length
            output_file = f'{self.output_dir}/{self.sample_name}.{x}.fa'
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
        # key: peptide value: list of ids
        fasta_df = self.create_index_file(unique_kmers)
        self.create_epitope_fastas(fasta_df)
