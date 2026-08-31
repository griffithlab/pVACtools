import sys
import argparse
import os
import shutil

from pvactools.lib.create_peptide_ordering_form import PvacfuseCreatePeptideOrderingForm

def define_parser():
    return PvacfuseCreatePeptideOrderingForm.parser('pvacfuse')

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    params = {
        'transcripts_fasta': args.transcripts_fasta,
        'flanking_sequence_length': args.flanking_sequence_length,
        'classI_aggregated_tsv': args.classI_aggregated_tsv,
        'classII_aggregated_tsv': args.classII_aggregated_tsv,
        'output_file_prefix': args.output_file_prefix,
        'sample_name': args.sample_name,
        'output_path': args.output_path,
        'aggregate_report_evaluation': args.aggregate_report_evaluation,
        'classI_IC50': args.classI_IC50,
        'classI_percent': args.classI_percent,
        'classII_IC50': args.classII_IC50,
        'classII_percent': args.classII_percent,
        'prob_pos': args.prob_pos,
    }
    creator = PvacfuseCreatePeptideOrderingForm(**params)
    creator.execute()

if __name__ == "__main__":
    main()
