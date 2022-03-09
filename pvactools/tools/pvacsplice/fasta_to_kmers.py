import os
import re
import pandas as pd
from pyfaidx import Fasta

class FastaToKmers():
    def __init__(self, **kwargs):
        self.tscript_fasta   = Fasta(kwargs['fasta'])
        self.output_dir      = kwargs['output_dir']
        self.epitope_lengths = [8] #kwargs['epitope_lengths']
        self.combined_df     = kwargs['combined_df']
        self.unique_kmers = {}


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
        # all kmers not in wt kmers
        final_kmers = {k: v for k, v in mut_dict.items() if k not in list(wt_dict.keys())}
        if len(final_kmers) == 0:
            print(f'{wt_name} does not produce any tumor-specific kmers.') # sys.exit()
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
        fasta_info = {k:':'.join(v) for k,v in self.unique_kmers.items()}
        # prepare data to convert to df
        modified_dict = {'kmer': fasta_info.keys(), 'indexes': fasta_info.values()}
        # convert to df
        self.index_df = pd.DataFrame.from_dict(modified_dict)
        # add length
        self.index_df['length'] = self.index_df['kmer'].str.len()
        # expand tscripts to one per line
        pd.options.mode.chained_assignment = None
        self.index_df['indexes'] = self.index_df.indexes.apply(lambda x: x.split(':'))
        self.index_df = self.index_df.explode('indexes')
        # add peptide_position
        self.index_df['peptide_position'] = self.index_df.indexes.apply(lambda x: x.split('.')[-1])
        self.index_df['indexes'] = self.index_df.indexes.apply(lambda x: '.'.join(x.split('.')[:-1]))
        # save to_csv()
        self.index_df.to_csv(f'{self.output_dir}/kmer_index.tsv' ,sep='\t', index=False)

    def create_epitope_fastas(self):
        # for only 1 length at a time
        for x in self.epitope_lengths:
            index_subset = self.index_df[self.index_df['length'] == x]
            # 1 file per kmer length
            file = f'{self.output_dir}/epitope_length_{x}.fa'
            # loop over rows in subset df
            for row in index_subset.itertuples():
                print(row)
                write_str = f'>{row.indexes}\n{row.kmer}\n'
                # don't duplicate entries
                if os.path.exists(file):
                    dup_content = re.search(write_str, open(file, 'r').read())
                    if dup_content == None:
                        with open(file, 'a') as f:
                            f.write(write_str)
                else:
                    with open(file, 'a') as f:
                        f.write(write_str)  

    def execute(self):
        self.loop_through_tscripts()
        self.create_index_file()
        self.create_epitope_fastas()
