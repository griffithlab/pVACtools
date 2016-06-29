import argparse
import sys
import re
import os
import csv

def coverage(ref, var):
    return ref + var

def vaf(ref, var):
    return (var / (coverage(ref, var)+0.00001)) * 100

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq coverage_filter')
    parser.add_argument('input_file', type=argparse.FileType('r'),
                        help="Output file from binding filter with read count info appended")
    parser.add_argument('output_file', type=argparse.FileType('w'),
                        help="Filtered file")
    parser.add_argument('--normal-cov', type=int,
                        help="Normal Coverage Cutoff. Sites above this cutoff will be considered. " +
                        "default 5",
                        default=5)
    parser.add_argument('--tdna-cov', type=int,
                        help="Tumor DNA Coverage Cutoff. Sites above this cutoff will be considered. " +
                        "default 10",
                        default=10)
    parser.add_argument('--trna-cov', type=int,
                        help="Tumor RNA Coverage Cutoff. Sites above this cutoff will be considered. " +
                        "default 10",
                        default=10)
    parser.add_argument('--normal-vaf', type=int,
                        help="Normal VAF Cutoff. Sites BELOW this cutoff in normal will be considered. " +
                        "default 2",
                        default=2)
    parser.add_argument('--tdna-vaf', type=int,
                        help="Tumor DNA VAF Cutoff. Sites above this cutoff will be considered. " +
                        "default 40",
                        default=40)
    parser.add_argument('--trna-vaf', type=int,
                        help="Tumor RNA VAF Cutoff. Sites above this cutoff will be considered. " +
                        "default 40",
                        default=40)
    parser.add_argument('--expn-val', type=int,
                        help="Gene Expression (FPKM) Cutoff. " +
                        "default 1",
                        default=1)
    args = parser.parse_args(args_input)

#### INPUT AND OUTPUT FILE FORMAT ##
#Chromosome
#Start
#Stop
#Reference
#Variant
#Transcript
#Ensembl Gene ID
#Variant Type
#Mutation
#Protein Position
#Gene Name
#HLA Allele
#Peptide Length
#Sub-peptide Position
#MT score
#WT score
#MT epitope seq
#WT epitope seq
#Fold Change
#Normal Ref Count
#Normal Var Count
#Tumor DNA Ref Count
#Tumor DNA Var Count
#Tumor RNA Ref Count
#Tumor RNA Var Count
#Gene Exp FPKM

    reader = csv.DictReader(args.input_file, delimiter='\t')
    writer = csv.DictWriter(
        args.output_file,
        reader.fieldnames,
        delimiter = '\t',
        lineterminator = '\n'
    )
    writer.writeheader()
    for entry in reader:
        if 'Normal Ref Count' in entry and 'Normal Var Count' in entry:
            ref = int(entry['Normal Ref Count'])
            var = int(entry['Normal Var Count'])
            ncov  = coverage(ref, var)
            if ncov < args.normal_cov:
                continue

            nvaf  = vaf(ref, var)
            if nvaf > args.normal_vaf:
                continue

        if 'Tumor DNA Ref Count' in entry and 'Tumor DNA Var Count' in entry:
            ref = int(entry['Tumor DNA Ref Count'])
            var = int(entry['Tumor DNA Var Count'])
            tdcov = coverage(ref, var)
            if tdcov < args.tdna_cov:
                continue

            tdvaf = vaf(ref, var)
            if tdvaf < args.tdna_vaf:
                continue

        if 'Tumor RNA Ref Count' in entry and 'Tumor RNA Var Count' in entry:
            ref = int(entry['Tumor RNA Ref Count'])
            var = int(entry['Tumor RNA Var Count'])
            trcov = coverage(ref, var)
            if trcov < args.trna_cov:
                continue

            trvaf = vaf(ref, var)
            if trvaf < args.trna_vaf:
                continue

        if 'Gene Exp FPKM' in entry:
            fpkm  = float(entry['Gene Exp FPKM'])
            if fpkm < args.expn_val:
                continue

        writer.writerow(entry)

    args.input_file.close()
    args.output_file.close()

if __name__ == "__main__":
    main()
