import argparse
import os
from subprocess import run, call, PIPE
import sys

current_directory = os.path.abspath(os.path.dirname(__file__))

def run_pvac_script(script, *args):
    pvac_cmd = "%s %s %s" % (
        sys.executable,
        os.path.join(current_directory, script),
        " ".join(str(arg) for arg in args)
    )
    result = call([pvac_cmd], shell=True)
    if result:
        print("Error while running:", script, args)
        exit(result)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("input",
                        help="Input TSV File with variants (please provide complete path)"
                        )
    parser.add_argument("sample_name",
                        help="Name of Sample; will be used as prefix for output files"
                        )
    parser.add_argument("netmhc_path",
                        help="Path to local NetMHC3.4 installation"
                        )
    parser.add_argument("allele",
                        help="Allele name to predict epitope prediction. Multiple alleles can be specified using a comma-separated list. For a list of available alleles, use: netMHC -A"
                        )
    parser.add_argument("epitope_length",
                        help="length of subpeptides(epitopes) to predict ; Multiple lengths can be specified using a comma-separated list. Typical epitope lengths vary between 8-11."
                        )
    parser.add_argument("output_dir",
                        help="Output directory for writing all result files"
                        )
    parser.add_argument("-l", "--variant-peptide-length",
                        type=int,
                        help="length of the peptide sequences in the input FASTA file; default 21",
                        default=21,
                        dest="peptide_sequence_length"
                        )
    parser.add_argument("-b","--binding-threshold",
                        type=int,
                        help="report only epitopes where the mutant allele has ic50 binding scores below this value ; default 500",
                        default=500,
                        dest="binding_threshold"
                        )
    parser.add_argument("-c", "--fold-change",
                        type=int,
                        help="Minimum fold change between mutant binding score and wild-type score. The default is 0, which filters no results, but 1 is often a sensible default (requiring that binding is better to the MT than WT)",
                        default=0,
                        dest="minimum_fold_change"
                        )

    args = parser.parse_args()
    args.epitope_length = [int(epl) for epl in args.epitope_length.split(',')]
    args.allele = [a for a in args.allele.split(',')]

    fasta_file = args.sample_name + "_" + str(args.peptide_sequence_length) + ".fa"
    fasta_key_file = args.sample_name + "_" + str(args.peptide_sequence_length) + ".key"

    print("Generating Variant Peptide FASTA File")
    run_pvac_script("generate_variant_sequences.py",
                    args.input,
                    args.peptide_sequence_length,
                    os.path.join(args.output_dir, fasta_file)
                    )
    print("Completed")

    print("Generating FASTA Key File")
    run_pvac_script("generate_fasta_key.py",
                    os.path.join(args.output_dir, fasta_file),
                    os.path.join(args.output_dir, fasta_key_file)
                    )
    print("Completed")

    netmhc_output = set(
        run([args.netmhc_path+" -A"],
            shell=True, stdout=PIPE
            ).stdout.decode().split('\n')
        )

    for epl in args.epitope_length:
        for a in args.allele:
            if a in netmhc_output:
                net_out = ".".join([args.sample_name, a, str(epl), "netmhc.xls"])
                print("Running NetMHC on Allele", a,
                      "and Epitope Length", epl
                      )
                netmhc_cmd = "%s -a %s -l %d %s -x %s" % (
                    args.netmhc_path, a, epl,
                    os.path.join(args.output_dir, fasta_file),
                    os.path.join(args.output_dir, net_out)
                )
                call([netmhc_cmd], shell=True)
            else:
                print("NetMHC allele not valid.  Please check using:")
                print(args.netmhc_path,"-A")
            print("Completed")

    fof = os.path.join(args.output_dir, args.sample_name+".fof")

    writer = open(fof, mode='w')

    for epl in args.epitope_length:
        for a in args.allele:
            net_out = ".".join([args.sample_name, a, str(epl), "netmhc.xls"])
            net_parsed = ".".join([args.sample_name, a, str(epl), "netmhc.parsed"])
            print("Parsing NetMHC Output")
            run_pvac_script("parse_output_netmhc.py",
                            os.path.join(args.output_dir, net_out),
                            os.path.join(args.output_dir, fasta_key_file),
                            os.path.join(args.output_dir, net_parsed)
                            )
            print("Completed")
            writer.write(os.path.join(args.output_dir, net_parsed))
            writer.write("\n")
    writer.close()

    filt_out = os.path.join(args.output_dir, args.sample_name+"_filtered.xls")
    print("Running Binding Filters")
    run_pvac_script("binding_filter.py",
                    args.input, fof, filt_out,
                    '-c', args.minimum_fold_change,
                    '-b', args.binding_threshold
                    )
    print("Completed")
    print("\n")
    print("Done: pVac-Seq.py has completed. File", filt_out,
          "contains list of binding-filtered putative neoantigens")
    print("We recommend appending coverage information and running CoverageFilters.py to filter based on sequencing coverage information")


if __name__ == '__main__':
    main()
