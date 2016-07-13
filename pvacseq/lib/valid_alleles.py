import sys
from pathlib import Path # if you haven't already done so
root = str(Path(__file__).resolve().parents[1])
sys.path.append(root)

import argparse
from lib import pvacseq_utils

def prediction_method_lookup(prediction_method):
    prediction_method_lookup_dict = pvacseq_utils.prediction_method_to_iedb_lookup_dict()
    return prediction_method_lookup_dict[prediction_method]

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq valid_alleles')
    parser.add_argument('-p', '--prediction_algorithm',
                        choices=pvacseq_utils.prediction_methods(),
                        help="Print available alleles for this epitope prediction algorithm",
                        required=False)
    args = parser.parse_args(args_input)
    if args.prediction_algorithm is None:
        print('\n'.join(sorted(pvacseq_utils.valid_allele_names())))
    else:
        iedb_method = prediction_method_lookup(args.prediction_algorithm)
        print("\n".join(sorted(pvacseq_utils.valid_allele_names_for_method(iedb_method))))

if __name__ == "__main__":
    main()
