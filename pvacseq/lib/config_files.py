import sys
import argparse
from collections import OrderedDict

def additional_input_file_list_options():
    return OrderedDict([
        ('gene_expn_file', 'genes.fpkm_tracking file from Cufflinks'),
        ('transcript_expn_file', 'isoforms.fpkm_tracking file from Cufflinks'),
    ])

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq config_files')
    parser.add_argument("config_file_type",
                        choices=['additional_input_file_list'],
                        help="The config file type to retrieve more information for",
                        )
    args = parser.parse_args(args_input)

    if args.config_file_type == 'additional_input_file_list':
        for item, description in additional_input_file_list_options().items():
            print("%s: <%s>" % (item, description))

if __name__ == "__main__":
    main()
