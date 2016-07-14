import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import argparse
import requests
import re
from lib import pvacseq_utils

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq call_iedb')
    parser.add_argument('input_file', type=argparse.FileType('r'),
                        help="Input FASTA file")
    parser.add_argument('output_file', type=argparse.FileType('w'),
                        help="Output file from iedb")
    parser.add_argument('method',
                        choices=pvacseq_utils.iedb_prediction_methods(),
                        help="The iedb analysis method to use")
    parser.add_argument('allele',
                        help="Allele for which to make prediction")
    parser.add_argument('epitope_length', type=int, choices=[8,9,10,11,12,13,14,15],
                        help="Length of subpeptides (epitopes) to predict")
    args = parser.parse_args(args_input)
    pvacseq_utils.check_alleles_valid([args.allele])
    pvacseq_utils.check_allele_valid_for_method(args.allele, args.method)
    pvacseq_utils.check_length_valid_for_allele_and_method(args.epitope_length, args.allele, args.method)

    data = {
        'sequence_text': args.input_file.read(),
        'method':        args.method,
        'allele':        args.allele,
        'length':        args.epitope_length,
    }

    response = requests.post('http://tools-api.iedb.org/tools_api/mhci/', data=data)
    if response.status_code == 500:
        #Retry once
        response = requests.post('http://tools-api.iedb.org/tools_api/mhci/', data=data)
    if response.status_code != 200:
        sys.exit("Error posting request to IEDB.\n%s" % response.text)
    args.output_file.write(response.text)

    args.input_file.close()
    args.output_file.close()

if __name__ == "__main__":
    main()
