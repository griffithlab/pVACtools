import sys
import os
import argparse
import re
from lib.prediction_class import *
import time
from subprocess import run, PIPE

def filter_response(response_text):
    lines = response_text.splitlines()
    remaining_lines = lines.copy()
    for line in lines:
        if line.startswith(b"allele"):
            return b"\n".join(remaining_lines)
        else:
            remaining_lines.pop(0)

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq call_iedb')
    parser.add_argument('input_file', type=argparse.FileType('r'),
                        help="Input FASTA file")
    parser.add_argument('output_file',
                        help="Output file from iedb")
    parser.add_argument('method',
                        choices=PredictionClass.prediction_methods(),
                        help="The iedb analysis method to use")
    parser.add_argument('allele',
                        help="Allele for which to make prediction")
    parser.add_argument('-l', '--epitope-length', type=int, choices=[8,9,10,11,12,13,14,15],
                        help="Length of subpeptides (epitopes) to predict")
    parser.add_argument(
        "-r", "--iedb-retries",type=int,
        default=5,
        help="Number of retries when making requests to the IEDB RESTful web interface. Must be less than or equal to 100."
             + "Default: 5"
    )
    parser.add_argument(
        "-e", "--iedb-executable-path",
        help="The executable path of the local IEDB install"
    )
    args = parser.parse_args(args_input)

    PredictionClass.check_alleles_valid([args.allele])
    prediction_class = getattr(sys.modules[__name__], args.method)
    prediction_class_object = prediction_class()
    prediction_class_object.check_allele_valid(args.allele)

    if isinstance(prediction_class_object, MHCI):
        prediction_class_object.check_length_valid_for_allele(args.epitope_length, args.allele)

    if args.epitope_length is None and prediction_class_object.needs_epitope_length:
        sys.exit("Epitope length is required for class I binding predictions")

    (response_text, output_mode) = prediction_class_object.predict(args.input_file, args.allele, args.epitope_length, args.iedb_executable_path, args.iedb_retries)

    tmp_output_file = args.output_file + '.tmp'
    tmp_output_filehandle = open(tmp_output_file, output_mode)
    tmp_output_filehandle.write(response_text)
    tmp_output_filehandle.close()
    os.replace(tmp_output_file, args.output_file)

    args.input_file.close()

if __name__ == "__main__":
    main()
