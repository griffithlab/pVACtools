import sys
import argparse
from collections import OrderedDict

def additional_input_file_list_options():
    return OrderedDict([
        ('gene_expn_file', 'genes.fpkm_tracking file from Cufflinks'),
        ('transcript_expn_file', 'isoforms.fpkm_tracking file from Cufflinks'),
        ('normal_snvs_coverage_file', 'bam-readcount output file for normal BAM and snvs'),
        ('normal_indels_coverage_file', 'bam-readcount output file for normal BAM and indels'),
        ('tdna_snvs_coverage_file', 'bam-readcount output file for tumor DNA BAM and snvs'),
        ('tdna_indels_coverage_file', 'bam-readcount output file for tumor DNA BAM and indels'),
        ('trna_snvs_coverage_file', 'bam-readcount output file for tumor RNA BAM and snvs'),
        ('trna_indels_coverage_file', 'bam-readcount output file for tumor RNA BAM and indels'),
    ])

def define_parser():
    parser = argparse.ArgumentParser('pvacseq config_files')
    parser.add_argument(
        "config_file_type",
        choices=['additional_input_file_list'],
        help="The config file type to retrieve more information for",
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.config_file_type == 'additional_input_file_list':
        for item, description in additional_input_file_list_options().items():
            print("%s: <%s>" % (item, description))

if __name__ == "__main__":
    main()
