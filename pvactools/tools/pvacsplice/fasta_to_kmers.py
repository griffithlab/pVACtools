import os
import sys
import re
import argparse
from pyfaidx import Fasta

class FastaToKmers():
    def __init__(self, tscript_fasta:str, output_dir:str, classI_lengths:str='8,9,10,11', classII_lengths:str='12,13,14,15,16,17,18'):
        self.tscript_fasta = Fasta(tscript_fasta)
        self.output_dir = output_dir
        self.kmer_lengths = classI_lengths +classII_lengths
        self.unique_kmers = {}
        

    def create_kmers(self, seq_name):
        kmer_dict = {}
        sequence = str(self.tscript_fasta[seq_name])
        # loop over sequence
        for i in range(len(sequence)):
            for x in self.kmer_lengths:
                # create kmer
                k = sequence[i:x+i]
                if len(k) == x:
                    # add or append to dict
                    if not k in kmer_dict.keys():
                        kmer_dict[k] = {'length':x, 'name': []}
                    kmer_dict[k]['name'].append(f'{seq_name}.{i}')
        return kmer_dict
    
    def save_kmer_dicts(self, wt_name, mt_name):
        wt_dict = self.create_kmers(wt_name)
        mut_dict = self.create_kmers(mt_name)
        final_kmers = {k:v for k, v in mut_dict.items() if k not in list(wt_dict.keys())}
        if len(final_kmers) == 0:
            print(f'{self.tscript_id} does not produce any tumor-specific kmers.')
        return final_kmers

    def loop_through_tscripts(self):
        fasta_keys = list(self.tscript_fasta.keys())
        unique_keys = set([x.split('.', 1)[1] for x in fasta_keys])
        for key in unique_keys:
            names = [x for x in fasta_keys if key in x]
            final_kmers = self.save_kmer_dicts(names[0], names[1])
            if not final_kmers:
                continue
            for k,v in final_kmers.items():
                if not k in self.unique_kmers.keys():
                    self.unique_kmers[k] = v
                else:
                    self.unique_kmers[k]['name'].extend(v['name'])

    def create_fasta_pvacbind(self):
        for x in self.kmer_lengths:
            # join transcripts by into string
            fasta_info = {k:':'.join(v['name']) for k,v in self.unique_kmers.items() if v['length'] == x}
            # 1 file per kmer length
            file = f'{self.output_dir}/pvacbind_{x}mers.fa'
            for k,v in fasta_info.items():
                write_str = f'>{v}\n{k}\n'
                if os.path.exists(file):
                    dup_content = re.search(write_str, open(file, 'r').read())
                    if dup_content == None:
                        with open(file, 'a') as f:
                            f.write(write_str)
                else:
                    with open(file, 'a') as f:
                        f.write(write_str)  


def define_parser():
    parser = argparse.ArgumentParser(
        description="Create pVACbind input fasta from predicted mutant transcript",
    )
    parser.add_argument('tscript_fasta', help='output fasta file with wt and mut aa seqs', type=str)
    parser.add_argument('output_dir', help='output path', type=str)
    parser.add_argument('classI_kmer_lengths', help='List of desired HLA classI neoantigen lengths', default='8,9,10,11', type=str)
    parser.add_argument('classII_kmer_lengths', help='List of desired HLA classII neoantigen lengths', default='12,13,14,15,16,17,18', type=str)
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    j = FastaToKmers(args.tscript_fasta, args.output_dir)
    j.loop_through_tscripts()
    j.create_fasta_pvacbind()
