import os
import sys
import yaml
import pkg_resources
from pvactools.lib.run_argument_parser import *

args_input = sys.argv[1:]
parser = PvacspliceRunArgumentParser().parser
args = parser.parse_args(args_input)

# Additional input: fasta_path - /Users/mrichters/Documents/git/pVACtools/tests/test_data/pvacsplice/inputs/all_sequences_chr1.fa
def print_log(log_dir='/Users/mrichters/Documents/current_testing/log', args_dict=vars(args)):
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, 'inputs.yaml')
    if os.path.exists(log_file):
        with open(log_file, 'r') as log_fh:
            past_inputs = yaml.load(log_fh, Loader=yaml.FullLoader)
            current_inputs = args_dict
            current_inputs['pvactools_version'] = pkg_resources.get_distribution("pvactools").version
            if past_inputs['pvactools_version'] != current_inputs['pvactools_version']:
                print('dif inputs!')
                sys.exit(
                    "Restart to be executed with a different pVACtools version:\n" +
                    "Past version: %s\n" % past_inputs['pvactools_version'] +
                    "Current version: %s" % current_inputs['pvactools_version']
                )
            for key in current_inputs.keys():
                if key == 'pvactools_version' or key == 'pvacseq_version':
                    continue
                if key not in past_inputs.keys() and current_inputs[key] is not None:
                    sys.exit(
                        "Restart inputs are different from past inputs: \n" +
                        "Additional input: %s - %s\n" % (key, current_inputs[key]) +
                        "Aborting."
                    )
                elif current_inputs[key] != past_inputs[key]:
                    print('key values dont match')
                    sys.exit(
                        "Restart inputs are different from past inputs: \n" +
                        "Past input: %s - %s\n" % (key, past_inputs[key]) +
                        "Current input: %s - %s\n" % (key, current_inputs[key]) +
                        "Aborting."
                    )
    else:
        with open(log_file, 'w') as log_fh:
            inputs = args_dict
            inputs['pvactools_version'] = pkg_resources.get_distribution("pvactools").version
            yaml.dump(inputs, log_fh, default_flow_style=False)

if __name__ == '__main__':
    print_log()