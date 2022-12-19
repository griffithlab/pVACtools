import os
import re
import pandas as pd
import numpy as np
from pyfaidx import Fasta

class FastaToKmers():
    def __init__(self, **kwargs):
        self.tscript_fasta   = Fasta(kwargs['fasta'])
        self.fasta_path    = kwargs['fasta']
        self.output_dir    = kwargs['output_dir']
        self.class_i_epitope_length  = kwargs['class_i_epitope_length']
        self.class_ii_epitope_length = kwargs['class_ii_epitope_length']
        self.class_i_hla   = kwargs['class_i_hla']
        self.class_ii_hla  = kwargs['class_ii_hla']
        self.sample_name   = kwargs['sample_name']
        # cumulative dict with all unique kmers from each transcript
        self.unique_kmers  = {}
        self.final_lengths = self.choose_final_lengths()

    def create_kmers(self, seq_name):
        kmer_dict = {}
        sequence = str(self.tscript_fasta[seq_name])
        # loop over entire sequence (doing this for WT and ALT) i == position in peptide
        for i in range(len(sequence)):
            # loop over lengths
            for x in self.final_lengths:
                final_seq_name = f'{seq_name};{i+1}'
                # grab kmer sequence
                k = sequence[i:x+i]
                # if kmer matches target len
                # if k seq doies not include X (any aa); continue
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
        final_kmers = {k: v for k, v in mut_dict.items() if k not in list(wt_dict.keys())}
        # for print statement
        junction_name = ' '.join(wt_name.split('.')[1:])
        if len(final_kmers) == 0:
            print(f'{junction_name} does not produce any tumor-specific kmers.')
            final_kmers = {}
        return final_kmers

    def loop_through_tscripts(self):
        # all fasta headers (WT and ALT)
        fasta_keys = list(self.tscript_fasta.keys())
        # take off WT. and ALT. prefixes from fasta_keys to de-duplicate
        unique_keys = set([x.split('.', 1)[1] for x in fasta_keys])
        for key in unique_keys:
            # selecting WT and ALT for each seq pair (bc calling fasta_keys not unique_keys)
            names = [x for x in fasta_keys if key in x]
            # get final mutated kmer list from save_kmer_dicts()
            final_kmers = self.create_unique_kmer_dict(names[0], names[1])
            if not final_kmers:
                continue
            # create master dict of unique kmers: index(es)
            for k,v in final_kmers.items():
                de_dup_v = v.replace('ALT.', '')
                if k not in self.unique_kmers.keys():
                    self.unique_kmers[k] = [de_dup_v]
                elif de_dup_v not in self.unique_kmers[k]:
                    self.unique_kmers[k].append(de_dup_v)

    def create_index_file(self):
        # joined indexes for each unique kmer - to format for df since lists are dif sizes
        fasta_info = {k:','.join(v) for k,v in self.unique_kmers.items()}
        # add header cols
        headers_dict = {'peptide': fasta_info.keys(), 'name': fasta_info.values()}
        # convert to df
        self.fasta_df = pd.DataFrame.from_dict(headers_dict)
        # save to sep dictionary
        index_df = self.fasta_df.copy()
        # convert string to list and explode so each index is on a sep line
        index_df['name'] = index_df['name'].str.split(',')
        index_df = index_df.explode('name')
        # now split each index into index, pos, len
        index_df[['junction_index', 'transcript_position']] = index_df['name'].str.split(';', expand=True)
        index_df['name'] = index_df['name'].str.replace(';', '.')
        # I don't think I need to save this file
        # create a tsv file
        # if os.path.exists(f'{self.output_dir}/kmer_index.tsv'):
        #    print('Kmer index already exists. Skipping.')
        # else:
        index_df.to_csv(f'{self.output_dir}/kmer_index.tsv', sep='\t', index=False)
        print('Kmer index file - complete')

    def choose_final_lengths(self):
        if not self.class_i_hla:
            lengths = self.class_ii_epitope_length
        elif not self.class_ii_hla:
            lengths = self.class_i_epitope_length
        return lengths

    def create_epitope_fastas(self):
        self.fasta_df['name'] = self.fasta_df['name'].str.replace(';', '.')
        # for only 1 length at a time
        for x in self.final_lengths:
            len_subset = self.fasta_df[self.fasta_df['peptide'].str.len() == x].sort_values(by=['name'])
            # 1 file per kmer length
            output_file = f'{self.output_dir}/{self.sample_name}.{x}.fa'
            # loop over rows in subset df
            for row in len_subset.itertuples():
                # fasta entry
                write_str = f'>{row.name}\n{row.peptide}\n'
                # don't duplicate entries
                if os.path.exists(output_file):
                    dup_content = re.search(row.peptide, open(output_file, "r").read())
                    if dup_content == None:
                        with open(output_file, "a") as f:
                            f.write(write_str)
                else:
                    with open(output_file, "a") as f:
                        f.write(write_str)
    
    def execute(self):
        self.loop_through_tscripts()
        self.create_index_file()
        self.create_epitope_fastas()
