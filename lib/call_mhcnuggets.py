import sys
import os
import tempfile
import argparse
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
stderr = sys.stderr
sys.stderr = open(os.devnull, 'w')
try:
    from mhcnuggets.src.predict import predict
except Exception as err:
    sys.stderr = stderr
    raise err
sys.stderr = stderr

def find_neoepitopes(sequence, length):
    epitopes = defaultdict(list)
    for i in range(0, len(sequence)-length+1):
        epitope = sequence[i:i+length]
        epitopes[epitope].append(i+1)
    return epitopes

def mhcnuggets_allele(allele, class_type):
    if class_type == 'I':
        return allele.replace('*', '')
    elif class_type == 'II':
        return "HLA-{}".format(allele).replace('*', '')

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('mhcnuggets', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file',
                        help="Input FASTA file")
    parser.add_argument('allele',
                        help="Allele for which to make prediction")
    parser.add_argument('epitope_length', type=int, choices=list(range(1,31)),
                        help="Length of subpeptides (epitopes) to predict")
    parser.add_argument('class_type', choices=['I','II'],
                        help="Class I or class II")
    parser.add_argument('output_file',
                        help="Output file from iedb")
    args = parser.parse_args(args_input)

    epitope_seq_nums = defaultdict(list)
    for record in SeqIO.parse(args.input_file, "fasta"):
        seq_num = record.id
        peptide = str(record.seq)
        epitopes = find_neoepitopes(peptide, args.epitope_length)
        for epitope, starts in epitopes.items():
            for start in starts:
                epitope_seq_nums[epitope].append((seq_num, start))

    tmp_file = tempfile.NamedTemporaryFile('w', delete=False)
    for epitope in epitope_seq_nums.keys():
        tmp_file.write("{}\n".format(epitope))
    tmp_file.close()

    tmp_output_file = tempfile.NamedTemporaryFile('r', delete=False)
    predict(args.class_type, tmp_file.name, mhcnuggets_allele(args.allele, args.class_type), output=tmp_output_file.name)
    os.unlink(tmp_file.name)
    tmp_output_file.close()
    df = pd.read_csv(tmp_output_file.name)
    os.unlink(tmp_output_file.name)
    processed_df = pd.DataFrame()
    for index, row in df.iterrows():
        seq_nums = epitope_seq_nums[row['peptide']]
        for seq_num, start in seq_nums:
            new_row = row.copy()
            new_row['seq_num'] = seq_num
            new_row['start'] = start
            new_row['allele'] = args.allele
            processed_df = processed_df.append(new_row)
    processed_df['start'] = pd.to_numeric(processed_df['start'], downcast='integer')
    processed_df = processed_df[['peptide', 'ic50', 'seq_num', 'start', 'allele']]
    processed_df.to_csv(args.output_file, index=False)

if __name__ == "__main__":
    main()
