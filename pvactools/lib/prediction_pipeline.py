import os

from pvactools.lib.call_predictors import CallPredictors
from pvactools.lib.combine_parsed_outputs import CombineParsedOutputs

class PredictionPipeline:
    def __init__(self, **kwargs):
        self.base_output_dir = kwargs['base_output_dir']
        self.mhc_class = kwargs['mhc_class']
        self.output_dir = os.path.join(self.base_output_dir, 'MHC_Class_{}'.format(self.mhc_class))
        self.sample_name = kwargs['sample_name']
        self.fasta_size = kwargs['fasta_size']
        self.prediction_algorithms = kwargs['prediction_algorithms']
        self.alleles = kwargs['alleles']
        self.epitope_lengths = kwargs['epitope_lengths']
        self.iedb_executable = kwargs['iedb_executable']
        self.iedb_retries = kwargs['iedb_retries']
        self.n_threads = kwargs['n_threads']
        self.use_normalized_percentiles = kwargs['use_normalized_percentiles']
        self.reference_scores_path = kwargs['reference_scores_path']
        self.additional_report_columns = kwargs['additional_report_columns']
        self.transcript_fasta = kwargs.pop('transcript_fasta', None)
        self.input_tsv_file = os.path.join(self.base_output_dir, f'{self.sample_name}.tsv')
        self.output_files = []

    def execute(self):
        print("Executing MHC Class {} predictions".format(self.mhc_class))

        os.makedirs(self.output_dir, exist_ok=True)

        for epitope_length in self.epitope_lengths:
            per_epitope_output_dir = os.path.join(self.output_dir, str(epitope_length))
            os.makedirs(per_epitope_output_dir, exist_ok=True)
            input_file = os.path.join(self.base_output_dir, f'{self.sample_name}.{epitope_length}.fa')
            if os.path.getsize(input_file) == 0:
                print("The intermediate FASTA file for epitope length {} is empty. No processable variants found.")
                continue

            parsed_output_files = []
            for allele in self.alleles:
                predictor_arguments = {
                    'input_file': input_file,
                    'sample_name': self.sample_name,
                    'fasta_size': self.fasta_size,
                    'allele': allele,
                    'epitope_length': epitope_length,
                    'prediction_algorithms': self.prediction_algorithms,
                    'iedb_executable_path': self.iedb_executable,
                    'iedb_retries': self.iedb_retries,
                    'n_threads': self.n_threads,
                    'output_dir': per_epitope_output_dir,
                }
                call_predictors = CallPredictors(**predictor_arguments)
                call_predictors.execute()

                if len(call_predictors.output_files) > 0:
                    parsed_file_path = os.path.join(call_predictors.tmp_dir, f"{self.sample_name}.{allele}.{epitope_length}.parsed.tsv")
                    if os.path.exists(parsed_file_path):
                        print(f"Parsed Output File for Allele {allele} and Epitope Length {epitope_length} already exists. Skipping")
                        parsed_output_files.append(parsed_file_path)
                        continue
                    print(f"Parsing prediction file for Allele {allele} and Epitope Length {epitope_length}")
                    parser_arguments = {
                        'prediction_files'          : call_predictors.output_files,
                        'tsv_file'                  : self.input_tsv_file,
                        'key_files'                 : call_predictors.output_key_files,
                        'output_file'               : parsed_file_path,
                        'use_normalized_percentiles': self.use_normalized_percentiles,
                        'reference_scores_path'     : self.reference_scores_path,
                        'sample_name'               : self.sample_name,
                    }
                    if self.additional_report_columns and 'sample_name' in self.additional_report_columns:
                        parser_arguments['add_sample_name_column'] = True
                    self.call_parser(parser_arguments)
                    print(f"Parsing prediction file for Allele {allele} and Epitope Length {epitope_length} - Completed")
                    parsed_output_files.append(parsed_file_path)
                else:
                    print(f"No predictions made for Allele {allele} and Epitope Length {epitope_length}")

            if len(parsed_output_files) > 0:
                print("Combining Parsed Prediction Files")
                output_file = os.path.join(per_epitope_output_dir, "{}.all_epitopes.tsv".format(self.sample_name))
                params = {
                    'input_files': parsed_output_files,
                    'output_file': output_file,
                }
                CombineParsedOutputs(**params).execute()
                print("Completed")

                if os.path.exists(output_file):
                    self.output_files.append(output_file)

                self.generate_fasta()
            else:
                print(f"No output files created for MHC Class {self.mhc_class}. Aborting")
