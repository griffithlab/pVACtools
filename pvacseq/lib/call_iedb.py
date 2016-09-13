import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import argparse
import requests
import re
import os
from lib.prediction_class import *

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
    args = parser.parse_args(args_input)

    PredictionClass.check_alleles_valid([args.allele])
    prediction_class = PredictionClass.prediction_class_for_iedb_prediction_method(args.method)
    prediction_class.check_allele_valid(args.allele)
    prediction_class_object = PredictionClass.prediction_class_for_iedb_prediction_method(args.method)

    if isinstance(prediction_class_object, MHCI):
        prediction_class.check_length_valid_for_allele(args.epitope_length, args.allele)

    if args.epitope_length is None and isinstance(prediction_class_object, MHCI):
        sys.exit("Epitope length is required for class I binding predictions")

    data = {
        'sequence_text': args.input_file.read(),
        'method':        args.method,
        'allele':        args.allele,
    }
    if args.epitope_length is not None:
        data['length'] = args.epitope_length

    url = prediction_class_object.url

    response = requests.post(url, data=data)
    if response.status_code == 500:
        #Retry once
        response = requests.post(url, data=data)
        if response.status_code == 500:
            #Retry a second time
            response = requests.post(url, data=data)
    if response.status_code != 200:
        sys.exit("Error posting request to IEDB.\n%s" % response.text)

    tmp_output_file = args.output_file + '.tmp'
    tmp_output_filehandle = open(tmp_output_file, 'w')
    tmp_output_filehandle.write(response.text)
    tmp_output_filehandle.close()
    os.replace(tmp_output_file, args.output_file)

    args.input_file.close()

if __name__ == "__main__":
    main()
