import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import argparse
import os
from subprocess import run, call, PIPE

try:
    from .. import lib
except ValueError:
    import lib 
from lib import pvacseq_utils

def prediction_method_lookup(prediction_method):
    prediction_method_lookup_dict = pvacseq_utils.prediction_method_to_iedb_lookup_dict()
    return prediction_method_lookup_dict[prediction_method]

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser("pvacseq run")

    parser.add_argument("input_file",
                        help="Input VCF with VEP annotations (please provide complete path)"
                        )
    parser.add_argument("sample_name",
                        help="Name of Sample; will be used as prefix for output files"
                        )
    parser.add_argument("allele",
                        help="Allele name to predict epitope prediction. Multiple alleles can be specified using a comma-separated list. For a list of available alleles, use: netMHC -A",
                        type=lambda s:[a for a in s.split(',')]
                        )
    parser.add_argument("epitope_length",
                        help="length of subpeptides(epitopes) to predict ; Multiple lengths can be specified using a comma-separated list. Typical epitope lengths vary between 8-11.",
                        type=lambda s:[int(epl) for epl in s.split(',')]
                        )
    parser.add_argument("prediction_algorithms",
                        choices=pvacseq_utils.prediction_methods(),
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
    parser.add_argument("-b","--binding-threshold",
                        type=int,
                        help="report only epitopes where the mutant allele has ic50 binding scores below this value ; default 500",
                        default=500)
    parser.add_argument("-c", "--minimum-fold-change",
                        type=int,
                        help="Minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default (requiring that binding is better to the MT than WT)",
                        default=0)

    args = parser.parse_args(args_input)

    print("Converting VCF to TSV")
    tsv_file = args.sample_name + '.tsv'
    lib.convert_vcf.main(
        [
            args.input_file,
            os.path.join(args.output_dir, tsv_file),
        ]
    )
    print("Completed")

    print("Generating Variant Peptide FASTA File")
    fasta_file = args.sample_name + "_" + str(args.peptide_sequence_length) + ".fa"
    fasta_key_file = args.sample_name + "_" + str(args.peptide_sequence_length) + ".key"
    lib.generate_fasta.main(
        [
            os.path.join(args.output_dir, tsv_file),
            str(args.peptide_sequence_length),
            os.path.join(args.output_dir, fasta_file)
        ]
    )
    print("Completed")

    print("Generating FASTA Key File")
    lib.generate_fasta_key.main(
        [
            os.path.join(args.output_dir, fasta_file),
            os.path.join(args.output_dir, fasta_key_file)
        ]
    )
    print("Completed")

#    netmhc_output = set(
#        run([args.netmhc_path+" -A"],
#            shell=True, stdout=PIPE
#            ).stdout.decode().split('\n')
#        )

    iedb_output_files = []
    for method in args.prediction_algorithms:
        iedb_method = prediction_method_lookup(method)
        for epl in args.epitope_length:
            for a in args.allele:
 #               if a in netmhc_output:
                iedb_out = os.path.join(args.output_dir, ".".join([args.sample_name, a, str(epl), iedb_method, "tsv"]))
                print("Running IEDB on Allele %s and Epitope Length %s with Method %s" % (a, epl, method))
                lib.call_iedb.main([
                    os.path.join(args.output_dir, fasta_file),
                    iedb_out,
                    iedb_method,
                    a,
                    str(epl),
                ])
                iedb_output_files.append(iedb_out)
#                else:
#                    print("Allele not valid.  Please check using:")
#                    print(args.netmhc_path,"-A")
                print("Completed")

    input_files = []

    for epl in args.epitope_length:
        for a in args.allele:
            net_parsed = ".".join([args.sample_name, a, str(epl), "parsed.tsv"])
            print("Parsing NetMHC Output")
            lib.parse_output.main(
                [
                    *iedb_output_files,
                    os.path.join(args.output_dir, tsv_file),
                    os.path.join(args.output_dir, fasta_key_file),
                    os.path.join(args.output_dir, net_parsed)
                ]
            )
            print("Completed")
            input_files.append(os.path.join(args.output_dir, net_parsed))

    filt_out = os.path.join(args.output_dir, args.sample_name+"_filtered.tsv")
    print("Running Binding Filters")
    lib.binding_filter.main(
        [
            *input_files,
            filt_out,
            '-c', str(args.minimum_fold_change),
            '-b', str(args.binding_threshold)
        ]
    )
    print("Completed")
    print("\n")
    print("Done: pvacseq has completed. File", filt_out,
          "contains list of binding-filtered putative neoantigens")
    print("We recommend appending coverage information and running CoverageFilters.py to filter based on sequencing coverage information")


if __name__ == '__main__':
    main()
