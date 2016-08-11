import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import argparse
import os
from subprocess import run, call, PIPE
import glob
import csv

try:
    from .. import lib
except ValueError:
    import lib
from lib.prediction_class import *
import shutil

def convert_vcf(args, output_dir):
    print("Converting VCF to TSV")
    tsv_file      = args.sample_name + '.tsv'
    tsv_file_path = os.path.join(output_dir, tsv_file)
    lib.convert_vcf.main([
        args.input_file,
        tsv_file_path,
    ])
    print("Completed")
    return tsv_file_path

def generate_fasta(args, tsv_file_path, output_dir):
    print("Generating Variant Peptide FASTA File")
    fasta_file      = args.sample_name + "_" + str(args.peptide_sequence_length) + ".fa"
    fasta_file_path = os.path.join(output_dir, fasta_file)
    lib.generate_fasta.main([
        tsv_file_path,
        str(args.peptide_sequence_length),
        fasta_file_path
    ])
    print("Completed")
    return fasta_file_path

def split_fasta_basename(args, tmp_dir):
    return os.path.join(tmp_dir, args.sample_name + "_" + str(args.peptide_sequence_length) + ".fa.split")

def split_fasta_file_and_create_key_files(args, fasta_file_path, tmp_dir):
    split_reader = open(fasta_file_path, mode='r')
    split_start = 1
    #Each fasta entry consists of two lines: header and sequence
    chunk_size  = args.fasta_size * 2
    chunks = []
    for chunk in split_file(split_reader, chunk_size):
        split_end = split_start + args.fasta_size - 1
        print("Splitting FASTA into smaller chunks - Entries %d-%d" % (split_start, split_end))
        split_fasta_file_path = "%s_%d-%d"%(split_fasta_basename(args, tmp_dir), split_start, split_end)
        if os.path.exists(split_fasta_file_path):
            print("Split FASTA file for Entries %d-%d already exists. Skipping." % (split_start, split_end))
            [entry for entry in chunk]
        else:
            split_writer = open(split_fasta_file_path, mode='w')
            split_writer.writelines(chunk)
            split_writer.close()
            print("Completed")
        print("Generating FASTA Key File - Entries %d-%d" % (split_start, split_end))
        split_fasta_key_file_path = split_fasta_file_path + '.key'
        if os.path.exists(split_fasta_key_file_path):
            print("Split FASTA Key File for Entries %d-%d already exists. Skipping." % (split_start, split_end))
        else:
            lib.generate_fasta_key.main([
                split_fasta_file_path,
                split_fasta_key_file_path,
            ])
            print("Completed")
        chunks.append("%d-%d" % (split_start, split_end))
        split_start += args.fasta_size
    split_reader.close()
    return chunks

def call_iedb_and_parse_outputs(args, chunks, tsv_file_path, tmp_dir):
    split_parsed_output_files = []
    for chunk in chunks:
        for a in args.allele:
            for epl in args.epitope_length:
                split_fasta_file_path = "%s_%s"%(split_fasta_basename(args, tmp_dir), chunk)
                split_iedb_output_files = []
                print("Processing entries for Allele %s and Epitope Length %s - Entries %s" % (a, epl, chunk))
                for method in args.prediction_algorithms:
                    prediction_class = globals()[method]
                    prediction = prediction_class()
                    iedb_method = prediction.iedb_prediction_method
                    valid_alleles = prediction.valid_allele_names()
                    if a not in valid_alleles:
                        print("Allele %s not valid for Method %s. Skipping." % (a, method))
                        continue
                    valid_lengths = prediction.valid_lengths_for_allele(a)
                    if epl not in valid_lengths:
                        print("Epitope Length %s is not valid for Method %s and Allele %s. Skipping." % (epl, method, a))
                        continue

                    split_iedb_out = os.path.join(tmp_dir, ".".join([args.sample_name, a, str(epl), iedb_method, "tsv_%s" % chunk]))
                    if os.path.exists(split_iedb_out):
                        print("IEDB file for Allele %s and Epitope Length %s with Method %s (Entries %s) already exists. Skipping." % (a, epl, method, chunk))
                        split_iedb_output_files.append(split_iedb_out)
                        continue
                    print("Running IEDB on Allele %s and Epitope Length %s with Method %s - Entries %s" % (a, epl, method, chunk))
                    lib.call_iedb.main([
                        split_fasta_file_path,
                        split_iedb_out,
                        iedb_method,
                        a,
                        str(epl),
                    ])
                    print("Completed")
                    split_iedb_output_files.append(split_iedb_out)

                split_parsed_file_path = os.path.join(tmp_dir, ".".join([args.sample_name, a, str(epl), "parsed", "tsv_%s" % chunk]))
                if os.path.exists(split_parsed_file_path):
                    print("Parsed Output File for Allele %s and Epitope Length %s (Entries %s) already exists. Skipping" % (a, epl, chunk))
                    split_parsed_output_files.append(split_parsed_file_path)
                    continue
                split_fasta_key_file_path = split_fasta_file_path + '.key'

                if len(split_iedb_output_files) > 0:
                    print("Parsing IEDB Output for Allele %s and Epitope Length %s - Entries %s" % (a, epl, chunk))
                    params = [
                        *split_iedb_output_files,
                        tsv_file_path,
                        split_fasta_key_file_path,
                        split_parsed_file_path,
                        '-m', args.top_score_metric,
                    ]
                    if args.top_result_per_mutation == True:
                        params.append('-t')
                    lib.parse_output.main(params)
                    print("Completed")
                    split_parsed_output_files.append(split_parsed_file_path)

    return split_parsed_output_files

def combined_parsed_outputs(args, split_parsed_output_files, output_dir):
    print("Combining Parsed IEDB Output Files")
    combined_parsed      = "%s.combined.parsed.tsv" % args.sample_name
    combined_parsed_path = os.path.join(output_dir, combined_parsed)
    lib.combine_parsed_outputs.main([
        *split_parsed_output_files,
        combined_parsed_path
    ])
    print("Completed")
    return combined_parsed_path

def binding_filter(args, combined_parsed_path, output_dir):
    filt_out_path = os.path.join(output_dir, args.sample_name+"_filtered.tsv")
    print("Running Binding Filters")
    lib.binding_filter.main(
        [
            combined_parsed_path,
            filt_out_path,
            '-c', str(args.minimum_fold_change),
            '-b', str(args.binding_threshold),
            '-m', str(args.top_score_metric),
        ]
    )
    print("Completed")
    return filt_out_path

def net_chop(args, input_path):
    output_path = os.path.join(args.output_dir, args.sample_name+"_filtered.chop.tsv")
    print("Submitting remaining epitopes to NetChop")
    lib.net_chop.main([
        input_path,
        output_path,
        '--method',
        args.net_chop_method,
        '--threshold',
        str(args.net_chop_threshold)
    ])
    print("Completed")
    return output_path

def netmhc_stab(args, input_path, output_dir):
    (filename, ext) = os.path.splitext(input_path)
    output_filepath = filename+'.stab'+ext
    print("Running NetMHCStabPan")
    lib.netmhc_stab.main(
        [
            input_path,
            output_filepath
        ]
    )
    print("Completed")
    return output_filepath

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser("pvacseq run")

    parser.add_argument("input_file",
                        help="Input VCF with VEP annotations (please provide complete path)"
                        )
    parser.add_argument("sample_name",
                        help="Name of Sample; will be used as prefix for output files"
                        )
    parser.add_argument("allele",
                        help="Allele name to predict epitope prediction. Multiple alleles can be specified using a comma-separated list. For a list of available alleles, use: pvacseq valid_alleles",
                        type=lambda s:[a for a in s.split(',')]
                        )
    parser.add_argument("epitope_length",
                        help="length of subpeptides(epitopes) to predict ; Multiple lengths can be specified using a comma-separated list. Typical epitope lengths vary between 8-11.",
                        type=lambda s:[int(epl) for epl in s.split(',')]
                        )
    parser.add_argument("prediction_algorithms",
                        choices=PredictionClass.prediction_methods(),
                        nargs="+",
                        help="The epitope prediction algorithms to use",
                        )
    parser.add_argument("output_dir",
                        help="Output directory for writing all result files"
                        )
    parser.add_argument("-l", "--peptide-sequence-length",
                        type=int,
                        help="length of the peptide sequences in the input FASTA file; default 21",
                        default=21)
    parser.add_argument('-t', '--top-result-per-mutation',
                        action='store_true',
                        help='Output top scoring candidate per allele-length per mutation. Default: False')
    parser.add_argument('-m', '--top-score-metric',
                        choices=['lowest', 'median'],
                        default='median',
                        help="The ic50 scoring metric to use when filtering epitopes. " +
                        "lowest: Best MT Score - lowest MT ic50 binding score of all chosen prediction methods. " +
                        "median: Median MT Score All Methods - median MT ic50 binding score of all chosen prediction methods. " +
                        "Default: median")
    parser.add_argument("-b","--binding-threshold",
                        type=int,
                        help="report only epitopes where the mutant allele has ic50 binding scores below this value ; default 500",
                        default=500)
    parser.add_argument("-c", "--minimum-fold-change",
                        type=int,
                        help="Minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default (requiring that binding is better to the MT than WT)",
                        default=0)
    parser.add_argument("-s", "--fasta-size",
                        type=int,
                        help="Number of fasta entries per IEDB request. For some resource-intensive prediction algorithms like Pickpocket and NetMHC it might be helpful to reduce this number. Needs to be an even number.",
                        default=200)
    parser.add_argument("-k", "--keep-tmp-files",
                        action='store_true',
                        help="Keep intermediate output files.",)
    parser.add_argument('--net-chop-method',
                        choices=lib.net_chop.methods,
                        help="NetChop prediction method to use (\"cterm\" for C term 3.0, \"20s\" for 20S 3.0).  Default: \"cterm\" (C term 3.0)",
                        default=None)
    parser.add_argument('--net-chop-threshold',
                        type=float,
                        help="NetChop prediction threshold.  Default: 0.5",
                        default=0.5)
    parser.add_argument(
        '--netmhc-stab',
        action='store_true',
        help="Run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes"
    )

    args = parser.parse_args(args_input)

    PredictionClass.check_alleles_valid(args.allele)

    if "." in args.sample_name:
        sys.exit("Sample name cannot contain '.'")

    if args.fasta_size%2 != 0:
        sys.exit("The fasta size needs to be an even number")

    output_dir = os.path.abspath(args.output_dir)

    tsv_file_path             = convert_vcf(args, output_dir)
    fasta_file_path           = generate_fasta(args, tsv_file_path, output_dir)

    if os.path.getsize(fasta_file_path) == 0:
        sys.exit("The fasta file is empty. Please check that the input VCF contains missense, inframe indel, or frameshift mutations.")

    tmp_dir = os.path.join(args.output_dir, 'tmp')
    os.makedirs(tmp_dir)
    chunks                    = split_fasta_file_and_create_key_files(args, fasta_file_path, tmp_dir)
    split_parsed_output_files = call_iedb_and_parse_outputs(args, chunks, tsv_file_path, tmp_dir)

    if len(split_parsed_output_files) == 0:
        sys.exit("No output files were created. Aborting.")

    combined_parsed_path      = combined_parsed_outputs(args, split_parsed_output_files, output_dir)
    final_path                = binding_filter(args, combined_parsed_path, output_dir)
    if args.net_chop_method:
        final_path = net_chop(
            args,
            final_path
        )

    if args.netmhc_stab:
        final_path = netmhc_stab(args, final_path, output_dir)

    print("\n")
    print("Done: pvacseq has completed. File", final_path,
          "contains list of binding-filtered putative neoantigens")
    print("We recommend appending coverage information and running `pvacseq coverage_filter` to filter based on sequencing coverage information")

    if not args.keep_tmp_files:
        shutil.rmtree(tmp_dir)

def split_file(reader, lines=400):
    from itertools import islice, chain
    tmp = next(reader)
    while tmp!="":
        yield chain([tmp], islice(reader, lines-1))
        try:
            tmp = next(reader)
        except StopIteration:
            return

if __name__ == '__main__':
    main()
