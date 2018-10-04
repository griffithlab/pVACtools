import sys
import os
import argparse
import re
from lib.prediction_class import *
from subprocess import run, PIPE

def setup_iedb_conda_env():
    env_check = run("conda env list | grep \"^pvactools_py27 \"", stdout=PIPE, shell=True)
    response = env_check.stdout.decode("utf-8")
    if response.count("\n") == 1:
        #environment with name "pvactools_py27" already exists; check that it really runs python2.7
        version_check = run("/bin/bash -c \"source activate pvactools_py27 && python -c \\\"import platform; print(platform.python_version())\\\"\"", stdout=PIPE, check=True, shell=True)
        if "2.7." not in version_check.stdout.decode("utf-8"):
            sys.exit('The existing conda environment "pvactools_py27" does not use python2.7. Please delete the existing environment.')
    elif response.count("\n") == 0:
        #environment with name "pvactools_py27" doesn't exist; create it
        run("conda create -n pvactools_py27 python=2.7 -y", check=True, shell=True)
    else:
        sys.exit("Something went wrong while checking the pvactools_py27 conda environment. `conda env list | grep \"^pvactools_py27 \"` returns more then one environment.")

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
    if output_mode == 'pandas':
        response_text.to_csv(tmp_output_file, index=False, sep="\t")
    else:
        tmp_output_filehandle = open(tmp_output_file, output_mode)
        tmp_output_filehandle.write(response_text)
        tmp_output_filehandle.close()
    os.replace(tmp_output_file, args.output_file)

    args.input_file.close()

if __name__ == "__main__":
    main()
