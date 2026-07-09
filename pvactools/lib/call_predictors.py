import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import pymp
import yaml
import logging

from pvactools.lib.prediction_class import *

class CallPredictors:
    def __init__(self, **kwargs):
        self.input_file = kwargs['input_file']
        self.sample_name = kwargs['sample_name']
        self.fasta_size = kwargs['fasta_size']
        self.allele = kwargs['allele']
        self.epitope_length = kwargs['epitope_length']
        self.prediction_algorithms = kwargs['prediction_algorithms'].copy()
        self.flurry_state = self.get_flurry_state()
        self.iedb_executable_path = kwargs['iedb_executable_path']
        self.iedb_retries = kwargs['iedb_retries']
        self.n_threads = kwargs['n_threads']
        self.output_dir = kwargs['output_dir']
        self.tmp_dir = os.path.join(self.output_dir, 'tmp')
        os.makedirs(self.tmp_dir, exist_ok=True)
        self.log_dir = os.path.join(self.output_dir, 'log')
        os.makedirs(self.log_dir, exist_ok=True),
        self.output_files = []
        self.output_key_files = []

    def get_flurry_state(self):
        if 'MHCflurry' in self.prediction_algorithms and 'MHCflurryEL' in self.prediction_algorithms:
            self.prediction_algorithms.remove('MHCflurryEL')
            return 'both'
        elif 'MHCflurry' in self.prediction_algorithms:
            return 'BA_only'
        elif 'MHCflurryEL' in self.prediction_algorithms:
            pred_idx = self.prediction_algorithms.index('MHCflurryEL')
            self.prediction_algorithms[pred_idx] = 'MHCflurry'
            return 'EL_only'
        else:
            return None

    def execute(self):
        split_fasta_files = self.split_input_file()

        argument_sets = []
        warning_messages = []
        for method in self.prediction_algorithms:
            prediction_class = globals()[method]
            prediction = prediction_class()
            valid_alleles = prediction.valid_allele_names()
            if self.allele not in valid_alleles:
                logging.info(f"Allele {self.allele} not valid for Method {method}. Skipping.")
                continue
            valid_lengths = prediction.valid_lengths_for_allele(self.allele)
            if self.epitope_length not in valid_lengths:
                logging.info(f"Epitope Length {self.epitope_length} is not valid for Method {method} and Allele {self.allele} Skipping.")
                continue
            for fasta_file in split_fasta_files:
                argument_sets.append([
                    fasta_file,
                    method,
                    self.allele
                ])
                _, chunk, _ = fasta_file.rsplit('.', 2)
                output_file = os.path.join(self.tmp_dir, f'{self.sample_name}.{method}.{self.allele}.{self.epitope_length}.{chunk}.tsv')
                self.output_files.append(output_file)

        with pymp.Parallel(self.n_threads) as p:
            for index in p.range(len(argument_sets)):
                arguments = argument_sets[index]
                filename = arguments[0]
                method = arguments[1]
                allele = arguments[2]
                p.print(f"Making binding predictions on Allele {allele} and Epitope Length {self.epitope_length} with Method {method} - File {filename}")

                prediction_class = getattr(sys.modules[__name__], method)
                prediction_class_object = prediction_class()

                (response_text, output_mode) = prediction_class_object.predict(filename, allele, self.epitope_length, self.iedb_executable_path, self.iedb_retries, tmp_dir=self.tmp_dir, log_dir=self.log_dir)

                _, chunk, _ = filename.rsplit('.', 2)
                output_file = os.path.join(self.tmp_dir, f'{self.sample_name}.{method}.{allele}.{self.epitope_length}.{chunk}.tsv')
                if output_mode == 'pandas':
                    response_text.to_csv(output_file, index=False, sep="\t")
                else:
                    with open(output_file, output_mode) as fh:
                        fh.write(response_text)

                p.print(f"Making binding predictions on Allele {allele} and Epitope Length {self.epitope_length} with Method {method} - File {filename} - Completed")

    def split_input_file(self):
        record_iterator = SeqIO.parse(self.input_file, "fasta")
        unique_records = defaultdict(list)
        for record in SeqIO.parse(self.input_file, "fasta"):
            unique_records[record.seq].append(record.id)

        file_index = 1
        current_records = []
        current_key_records = {}
        output_files = []
        for i, (key, value) in enumerate(unique_records.items()):
            start = (file_index - 1) * self.fasta_size + 1
            stop = file_index * self.fasta_size
            current_records.append(SeqRecord(key, id=str(i-start+2), description=""))
            current_key_records[i-start+2] = value

            # When we hit the specified number, write the file
            if len(current_records) == self.fasta_size:
                output_file = os.path.join(self.tmp_dir, f'{self.sample_name}.{start}-{stop}.fa')
                SeqIO.write(current_records, output_file, "fasta")

                output_key_file = os.path.join(self.tmp_dir, f'{self.sample_name}.{start}-{stop}.key')
                with open(output_key_file, 'w') as fh:
                    yaml.dump(current_key_records, fh, default_flow_style=False)

                # Reset for the next batch
                current_records = []
                current_key_records = {}
                file_index += 1
                output_files.append(output_file)
                self.output_key_files.append(output_key_file)

        if current_records:
            output_file = os.path.join(self.tmp_dir, f'{self.sample_name}.{start}-{stop}.fa')
            SeqIO.write(current_records, output_file, "fasta")

            output_key_file = os.path.join(self.tmp_dir, f'{self.sample_name}.{start}-{stop}.key')
            with open(output_key_file, 'w') as fh:
                yaml.dump(current_key_records, fh, default_flow_style=False)

            output_files.append(output_file)
            self.output_key_files.append(output_key_file)

        return output_files
