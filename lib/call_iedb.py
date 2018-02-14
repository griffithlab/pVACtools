import sys
import os
import argparse
import requests
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
                        choices=PredictionClass.iedb_prediction_methods(),
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
    prediction_class = PredictionClass.prediction_class_for_iedb_prediction_method(args.method)
    prediction_class.check_allele_valid(args.allele)
    prediction_class_object = PredictionClass.prediction_class_for_iedb_prediction_method(args.method)

    if isinstance(prediction_class_object, MHCI):
        prediction_class.check_length_valid_for_allele(args.epitope_length, args.allele)

    if args.epitope_length is None and prediction_class_object.needs_epitope_length:
        sys.exit("Epitope length is required for class I binding predictions")

    if args.iedb_executable_path is not None:
        response = run(prediction_class_object.iedb_executable_params(args), stdout=PIPE, check=True)
        response_text = filter_response(response.stdout)
        output_mode = 'wb'
    else:
        data = {
            'sequence_text': args.input_file.read(),
            'method':        args.method,
            'allele':        args.allele.replace('-DPB', '/DPB'),
            'user_tool':     'pVac-seq',
        }
        if args.epitope_length is not None:
            data['length'] = args.epitope_length

        url = prediction_class_object.url

        response = requests.post(url, data=data)
        retries = 0
        while response.status_code == 500 and retries < args.iedb_retries:
            time.sleep(60 * retries)
            response = requests.post(url, data=data)
            print("IEDB: Retry %s of %s" % (retries, args.iedb_retries))
            retries += 1

        if response.status_code != 200:
            sys.exit("Error posting request to IEDB.\n%s" % response.text)
        response_text = response.text
        output_mode = 'w'

    tmp_output_file = args.output_file + '.tmp'
    tmp_output_filehandle = open(tmp_output_file, output_mode)
    tmp_output_filehandle.write(response_text)
    tmp_output_filehandle.close()
    os.replace(tmp_output_file, args.output_file)

    args.input_file.close()

if __name__ == "__main__":
    main()
