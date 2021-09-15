import sys
import os
import argparse
import re
from subprocess import run, PIPE

from pvactools.lib.prediction_class import *

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq call_iedb', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file',
                        help="Input FASTA file")
    parser.add_argument('output_file',
                        help="Output file from iedb")
    parser.add_argument('method',
                        choices=PredictionClass.prediction_methods(),
                        help="The iedb analysis method to use")
    parser.add_argument('allele',
                        help="Allele for which to make prediction")
    parser.add_argument('-l', '--epitope-length', type=int, choices=[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30],
                        help="Length of subpeptides (epitopes) to predict")
    parser.add_argument(
        "-r", "--iedb-retries",type=int,
        default=5,
        help="Number of retries when making requests to the IEDB RESTful web interface. Must be less than or equal to 100."
    )
    parser.add_argument(
        "-e", "--iedb-executable-path",
        help="The executable path of the local IEDB install"
    )
    args = parser.parse_args(args_input)

    prediction_class = getattr(sys.modules[__name__], args.method)
    prediction_class_object = prediction_class()

    try:
        (response_text, output_mode) = prediction_class_object.predict(args.input_file, args.allele, args.epitope_length, args.iedb_executable_path, args.iedb_retries)
    except Exception as err:
        if str(err) == 'len(peptide_list) != len(scores)':
            (response_text, output_mode) = prediction_class_object.predict(args.input_file, args.allele, args.epitope_length, args.iedb_executable_path, args.iedb_retries)
        else:
            raise err

    tmp_output_file = args.output_file + '.tmp'
    if output_mode == 'pandas':
        response_text.to_csv(tmp_output_file, index=False, sep="\t")
    else:
        tmp_output_filehandle = open(tmp_output_file, output_mode)
        tmp_output_filehandle.write(response_text)
        tmp_output_filehandle.close()
    os.replace(tmp_output_file, args.output_file)

if __name__ == "__main__":
    main()
