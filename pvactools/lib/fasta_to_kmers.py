import os
import re
import pandas as pd
from pyfaidx import Fasta

class FastaToKmers():
    def __init__(self, **kwargs):
        self.tscript_fasta   = Fasta(kwargs['fasta'])
        self.fasta_path      = kwargs['fasta']
        self.output_dir      = kwargs['output_dir']
        self.epitope_lengths = kwargs['epitope_lengths']
        self.combined_df     = kwargs['combined_df']
        self.sample_name     = kwargs['sample_name']
        self.unique_kmers    = {}


    def create_kmers(self, seq_name):
        kmer_dict = {}
        sequence = str(self.tscript_fasta[seq_name])
        # loop over entire sequence (doing this for WT and ALT) i == position in peptide
        for i in range(len(sequence)):
            final_seq_name = f'{seq_name}.{i+1}'
            # loop over lengths
            for x in self.epitope_lengths:
                # grab kmer sequence
                k = sequence[i:x+i]
                # if kmer matches target len
                if len(k) == x:
                    # add entry to dictionary
                    kmer_dict[k] = final_seq_name
        return kmer_dict

    def save_kmer_dicts(self, wt_name, mt_name):
        wt_dict = self.create_kmers(wt_name)
        mut_dict = self.create_kmers(mt_name)
        junction_name = ' '.join(wt_name.split('.')[1:])
        # all kmers not in wt kmers
        final_kmers = {k: v for k, v in mut_dict.items() if k not in list(wt_dict.keys())}
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
            final_kmers = self.save_kmer_dicts(names[0], names[1])
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
        # joined indexes for each unique kmer
        fasta_info = {k:','.join(v) for k,v in self.unique_kmers.items()}
        # prepare data to convert to df
        modified_dict = {'kmer': fasta_info.keys(), 'indexes': fasta_info.values()}
        # convert to df
        # save this df to write to fasta file
        self.fasta_df = pd.DataFrame.from_dict(modified_dict)
        # add length
        self.fasta_df['length'] = self.fasta_df['kmer'].str.len()
        # expand tscripts to one per line in separate df and save to tsv
        index_df = self.fasta_df.copy()
        index_df['indexes'] = index_df.indexes.apply(lambda x: x.split(','))
        index_df = index_df.explode('indexes')
        if os.path.exists(f'{self.output_dir}/kmer_index.tsv'):
            print('Kmer index already exists. Skipping.')
        else:
            index_df.to_csv(f'{self.output_dir}/kmer_index.tsv' ,sep='\t', index=False)
            print('Kmer index file - complete')

    def create_epitope_fastas(self):
        # for only 1 length at a time
        for x in self.epitope_lengths:
            len_subset = self.fasta_df[self.fasta_df['length'] == x].sort_values(by=['indexes'])
            # 1 file per kmer length
            output_file = f'{self.output_dir}/{self.sample_name}.{x}.fa'
            # loop over rows in subset df
            for row in len_subset.itertuples():
                # fasta entry
                write_str = f'>{row.indexes}\n{row.kmer}\n'
                # don't duplicate entries
                if os.path.exists(output_file):
                    dup_content = re.search(write_str, open(output_file, "r").read())
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
