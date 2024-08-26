import os
import sys
import yaml
import pkg_resources
from pvactools.lib.run_argument_parser import *


def print_log(log_dir, args_dict, output_file_prefix):
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f'{output_file_prefix}.yml')
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
