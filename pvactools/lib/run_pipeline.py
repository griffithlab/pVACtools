import sys
import os
import platform
import logging

from pvactools.lib.prediction_class import NetMHCIIVersion
from pvactools.lib.print_log import *
from pvactools.lib.run_utils import *
from pvactools.lib.prediction_class_utils import *

class RunPipeline:
    def __init__(self, **kwargs):
        for (k,v) in kwargs.items():
            setattr(self, k, v)

        self.original_parameters = kwargs

        if self.iedb_retries > 100:
            raise Exception("The number of IEDB retries must be less than or equal to 100")

        if self.n_threads > 1 and platform.system() == "Darwin":
            raise Exception("Multithreading is not supported on MacOS")

        if (self.netmhciipan_version == '4.0' and self.iedb_install_directory is not None):
            raise Exception("Standalone IEDB does not support version 4.0")
        NetMHCIIVersion.netmhciipan_version = self.netmhciipan_version

        (self.class_i_prediction_algorithms, self.class_ii_prediction_algorithms) = split_algorithms(self.prediction_algorithms)
        alleles = combine_class_ii_alleles(self.allele)
        (self.class_i_alleles, self.class_ii_alleles, self.species) = split_alleles(alleles)

        if self.iedb_install_directory:
            self.iedb_mhc_i_executable = os.path.join(self.iedb_install_directory, 'mhc_i', 'src', 'predict_binding.py')
            if not os.path.exists(iedb_mhc_i_executable):
                raise Exception("IEDB MHC I executable path doesn't exist %s" % self.iedb_mhc_i_executable)
            self.iedb_mhc_ii_executable = os.path.join(self.iedb_install_directory, 'mhc_ii', 'mhc_II_binding.py')
            if not os.path.exists(iedb_mhc_ii_executable):
                raise Exception("IEDB MHC II executable path doesn't exist %s" % self.iedb_mhc_ii_executable)
        else:
            self.iedb_mhc_i_executable = None
            self.iedb_mhc_ii_executable = None

        if self.use_normalized_percentiles and self.species != 'human':
            logging.info("Normalized percentiles are only available for human alleles. Option will be ignored.")
            self.use_normalized_percentiles = False

        self.check_tumor_purity_argument()
        self.check_downstream_sequence_length_argument()
        self.extra_argument_checks()

        self.base_output_dir = os.path.abspath(self.output_dir)
        os.makedirs(self.base_output_dir, exist_ok=True)

    def check_tumor_purity_argument(self):
        if self.tumor_purity is not None:
            if self.tumor_purity > 1:
                raise Exception("--tumor-purity must be a float between 0 and 1. Value too large: {}".format(args.tumor_purity))
            elif self.tumor_purity < 0:
                raise Exception("--tumor-purity must be a float between 0 and 1. Value too small: {}".format(args.tumor_purity))

    def check_downstream_sequence_length_argument(self):
        if self.downstream_sequence_length == 'full':
            self.downstream_sequence_length = None
        elif self.downstream_sequence_length.isdigit():
            self.downstream_sequence_length = int(self.downstream_sequence_length)
        else:
            raise Exception("The downstream sequence length needs to be a positive integer or 'full'")

    def extra_argument_checks(self):
        pass

    def execute(self):
        print_log(os.path.join(self.base_output_dir, 'log'), self.original_parameters, 'inputs')
        self.call_input_to_kmer_pipeline()
        self.call_prediction_pipeline()
        self.call_ml_predictor()
        self.extra_steps()
        change_permissions_recursive(self.base_output_dir, 0o755, 0o644)

    def call_prediction_pipeline(self):
        all_params = {
            'I': {
                'iedb_executable': self.iedb_mhc_i_executable,
                'prediction_algorithms': self.class_i_prediction_algorithms,
                'alleles': self.class_i_alleles,
                'epitope_lengths': self.class_i_epitope_length,
                'netmhc_stab': self.netmhc_stab,
                'use_normalized_percentiles': self.use_normalized_percentiles,
                'reference_scores_path': self.reference_scores_path
            },
            'II': {
                'iedb_executable': self.iedb_mhc_ii_executable,
                'prediction_algorithms': self.class_ii_prediction_algorithms,
                'alleles': self.class_ii_alleles,
                'epitope_lengths': self.class_ii_epitope_length,
                'netmhc_stab': False,
                'use_normalized_percentiles': False,
                'reference_scores_path': self.reference_scores_path
            }
        }

        for (mhc_class, params) in all_params.items():
            if len(params['prediction_algorithms']) > 0 and len(params['alleles']) > 0:
                params['base_output_dir'] = self.base_output_dir
                params['mhc_class'] = mhc_class
                params['sample_name'] = self.sample_name
                params['fasta_size'] = self.fasta_size
                params['iedb_retries'] = self.iedb_retries
                params['n_threads'] = self.n_threads
                params['additional_report_columns'] = self.additional_report_columns
                params['transcript_fasta'] = self.transcript_fasta

                self.run_prediction_pipeline(params)

                if len(self.predictor.output_files) > 0:
                    all_epitopes_file = self.create_per_class_report(mhc_class)

                    post_processing_params = self.original_parameters.copy()
                    if self.run_reference_proteome_similarity:
                        post_processing_params['fasta'] = self.generate_flanked_protein_fasta(7)
                    # generate net_chop fasta to output dir if specified
                    if self.net_chop_method:
                        post_processing_params['net_chop_fasta'] = self.generate_flanked_protein_fasta(10)
                    filtered_file = os.path.join(self.predictor.output_dir, f"{self.sample_name}.MHC_{mhc_class}.filtered.tsv")
                    post_processing_params["filename_addition"] = f"MHC_{mhc_class}"
                    post_processing_params["species"] = self.species
                    post_processing_params["netmhc_stab"] = params["netmhc_stab"]
                    self.call_post_processor(all_epitopes_file, filtered_file, post_processing_params)
                else:
                    logging.info("\nNo processable variants found. Aborting.\n")
            elif len(params["prediction_algorithms"]) == 0:
                logging.info(f"No MHC class {mhc_class} prediction algorithms chosen. Skipping MHC class {mhc_class} predictions.")
            elif len(params["alleles"]) == 0:
                logging.info(f"No MHC class {mhc_class} alleles chosen. Skipping MHC class {mhc_class} predictions.")

    def generate_flanked_protein_fasta(self, flanking_sequence_length):
        fasta_file = os.path.join(self.predictor.output_dir, f"{self.sample_name}.{flanking_sequence_length}.fasta")
        if not os.path.exists(fasta_file):
            fasta_params = {
                'fasta_file_path': self.transcript_fasta,
                'trimmed_fasta_file_path': fasta_file,
                'flanking_sequence_length': flanking_sequence_length,
                'mutant_only': False,
            }
            self.call_generate_protein_fasta(fasta_params)
        return fasta_file

    def create_per_class_report(self, mhc_class):
        for file_name in self.predictor.output_files:
            if not os.path.exists(file_name):
                logging.info(f"File {self.file_name} doesn't exist. Aborting.")
                return
        all_epitopes_file = os.path.join(self.predictor.output_dir, f"{self.sample_name}.MHC_{mhc_class}.all_epitopes.tsv")
        combine_reports(self.predictor.output_files, all_epitopes_file)
        return all_epitopes_file

    def extra_steps(self):
        pass
