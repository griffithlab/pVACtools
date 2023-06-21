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
    parser.add_argument('--tmp-dir',
                        help="Location to write tmp files to")
    args = parser.parse_args(args_input)

    epitope_seq_nums = defaultdict(list)
    for record in SeqIO.parse(args.input_file, "fasta"):
        seq_num = record.id
        peptide = str(record.seq)
        epitopes = find_neoepitopes(peptide, args.epitope_length)
        for epitope, starts in epitopes.items():
            for start in starts:
                epitope_seq_nums[epitope].append((seq_num, start))

    tmp_file = tempfile.NamedTemporaryFile('w', dir=args.tmp_dir, delete=False)
    for epitope in epitope_seq_nums.keys():
        tmp_file.write("{}\n".format(epitope))
    tmp_file.close()

    tmp_output_file = tempfile.NamedTemporaryFile('w', dir=args.tmp_dir, delete=False)
    predict(args.class_type, tmp_file.name, mhcnuggets_allele(args.allele, args.class_type), output=tmp_output_file.name, rank_output=True)
    os.unlink(tmp_file.name)
    tmp_output_file.close()
    rank_output_file_name = "{}_ranks".format(tmp_output_file.name)
    os.unlink(tmp_output_file.name)

    df = pd.read_csv(rank_output_file_name)
    processed_df = pd.DataFrame()
    for index, row in df.iterrows():
        seq_nums = epitope_seq_nums[row['peptide']]
        for seq_num, start in seq_nums:
            new_row = row.copy().to_dict()
            new_row['seq_num'] = seq_num
            new_row['start'] = start
            new_row['allele'] = args.allele
            new_row['percentile'] = float(row['human_proteome_rank']) * 100
            processed_df = pd.concat([processed_df, pd.DataFrame.from_records(new_row, index=[0])], ignore_index=True)
    processed_df['start'] = pd.to_numeric(processed_df['start'], downcast='integer')
    processed_df = processed_df[['peptide', 'ic50', 'percentile', 'seq_num', 'start', 'allele']]
    processed_df.to_csv(args.output_file, index=False)
    os.unlink(rank_output_file_name)

if __name__ == "__main__":
    main()
