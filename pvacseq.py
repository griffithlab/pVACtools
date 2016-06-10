import argparse
import sys
from subprocess import call
import os

def pvac_path_for_script(script):
    pVac_directory = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(pVac_directory, "pvacseq", script)

def run_subcommand(script, args):
    script_cmd = "%s %s %s" % (
        sys.executable,
        pvac_path_for_script(script),
        " ".join(args)
    )
    call([script_cmd], shell=True)

def main():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()

    #add subcommands
    variant_sequences_parser = subparsers.add_parser("generate_fasta",
                                                     help="Generates a variant peptide FASTA file from the TSV input file",
                                                     add_help=False)
    variant_sequences_parser.set_defaults(script="generate_fasta.py")

    binding_filter_parser = subparsers.add_parser("binding_filter",
                                                  help="Filters variants processed by NetMHC",
                                                  add_help=False)
    binding_filter_parser.set_defaults(script="binding_filter.py")

    fasta_key_parser = subparsers.add_parser("generate_fasta_key",
                                             help="Generates a FASTA key file",
                                             add_help=False)
    fasta_key_parser.set_defaults(script="generate_fasta_key.py")

    parse_netmhc_parser = subparsers.add_parser("parse_output",
                                                help="Parses output from NetMHC",
                                                add_help=False)
    parse_netmhc_parser.set_defaults(script="parse_output.py")

    run_main_program_parser = subparsers.add_parser("run",
                                                    help="Runs the pVAC-Seq pipeline",
                                                    add_help=False)
    run_main_program_parser.set_defaults(script="main.py")

    args = parser.parse_known_args()
    try:
        run_subcommand(args[0].script, args[1])
    except AttributeError as e:
        parser.print_help()
        print("Error: No command specified")
        exit(-1)


if __name__ == '__main__':
    main()
