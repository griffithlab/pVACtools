from abc import ABCMeta, abstractmethod
import os
import csv
import sys
import inspect
from mhcflurry import Class1AffinityPredictor
from mhcnuggets.src.predict import predict
import requests
import re
import pandas as pd
import time
from subprocess import run, PIPE
import tempfile
from collections import defaultdict

class IEDB(metaclass=ABCMeta):
    @classmethod
    def iedb_prediction_methods(cls):
        return [prediction_class().iedb_prediction_method for prediction_class in cls.prediction_classes()]

    @abstractmethod
    def parse_iedb_allele_file(self):
        pass

    @abstractmethod
    def iedb_executable_params(self, args):
        pass

    @property
    @abstractmethod
    def iedb_prediction_method(self):
        pass

    @property
    @abstractmethod
    def url(self):
        pass

    @classmethod
    def filter_response(cls, response_text):
        lines = response_text.splitlines()
        remaining_lines = lines.copy()
        for line in lines:
            if line.startswith(b"allele"):
                return b"\n".join(remaining_lines)
            else:
                remaining_lines.pop(0)

    def predict(self, input_file, allele, epitope_length, iedb_executable_path, iedb_retries):
        if iedb_executable_path is not None:
            arguments = self.iedb_executable_params(iedb_executable_path, self.iedb_prediction_method, allele, input_file, epitope_length)
            response = run("/bin/bash -c \"source activate pvactools_py27; python {}\"".format(arguments), stdout=PIPE, check=True, shell=True)
            response_text = self.filter_response(response.stdout)
            return (response_text, 'wb')
        else:
            data = {
                'sequence_text': input_file.read(),
                'method':        self.iedb_prediction_method,
                'allele':        allele.replace('-DPB', '/DPB').replace('-DQB', '/DQB'),
                'user_tool':     'pVac-seq',
            }
            if epitope_length is not None:
                data['length'] = epitope_length

            response = requests.post(self.url, data=data)
            retries = 0
            while response.status_code == 500 and retries < iedb_retries:
                time.sleep(60 * retries)
                response = requests.post(self.url, data=data)
                print("IEDB: Retry %s of %s" % (retries, iedb_retries))
                retries += 1

            if response.status_code != 200:
                sys.exit("Error posting request to IEDB.\n%s" % response.text)
            response_text = response.text
            output_mode = 'w'
            return (response_text, 'w')

class PredictionClass(metaclass=ABCMeta):
    valid_allele_names_dict = {}
    allele_cutoff_dict = {}

    @classmethod
    def prediction_classes(cls):
        prediction_classes = []
        if not inspect.isabstract(cls):
            prediction_classes.append(cls)
        for subclass in cls.__subclasses__():
            prediction_classes.extend(subclass.prediction_classes())
        return prediction_classes

    @classmethod
    def prediction_methods(cls):
        return sorted([prediction_class.__name__ for prediction_class in cls.prediction_classes()])

    @classmethod
    def prediction_class_for_iedb_prediction_method(cls, method):
        prediction_classes = cls.prediction_classes()
        for prediction_class in prediction_classes:
            prediction_class_object = prediction_class()
            if ( issubclass(prediction_class_object.__class__, IEDBMHCI) or issubclass(prediction_class_object.__class__, IEDBMHCII) ) and prediction_class_object.iedb_prediction_method == method:
                return prediction_class_object
        module = getattr(sys.modules[__name__], method)
        return module()

    @classmethod
    def prediction_class_name_for_iedb_prediction_method(cls, method):
        return cls.prediction_class_for_iedb_prediction_method(method).__class__.__name__

    @classmethod
    def allele_info(cls, prediction_algorithms, name_filter):
        alleles = defaultdict(list)
        if prediction_algorithms is None:
            prediction_classes = cls.prediction_classes()
        else:
            prediction_classes = map(lambda a: globals()[a], prediction_algorithms.split(','))
        for prediction_class in prediction_classes:
            for allele in prediction_class().valid_allele_names():
                if name_filter is not None:
                    if name_filter.lower() in allele.lower():
                        alleles[allele].append(prediction_class.__name__)
                else:
                    alleles[allele].append(prediction_class.__name__)
        info = []
        for allele, prediction_algorithms in alleles.items():
            info.append({
                'name': allele,
                'prediction_algorithms': prediction_algorithms,
            })
        return info

    @classmethod
    def all_valid_allele_names(cls):
        valid_alleles = set()
        for prediction_class in cls.prediction_classes():
            valid_alleles.update(prediction_class().valid_allele_names())
        return list(valid_alleles)

    @classmethod
    def check_alleles_valid(cls, alleles):
        valid_alleles = cls.all_valid_allele_names()
        for allele in alleles:
            if allele not in valid_alleles:
                sys.exit("Allele %s not valid. Run `pvacseq valid_alleles` for a list of valid allele names." % allele)

    @classmethod
    def parse_allele_cutoff_file(cls):
        base_dir                = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        iedb_alleles_dir        = os.path.join(base_dir, 'tools', 'pvacseq', 'iedb_alleles')
        allele_cutoff_file_name = os.path.join(iedb_alleles_dir, "cutoffs.csv")
        cutoffs = {}
        with open(allele_cutoff_file_name) as allele_cutoff_file:
            csv_reader = csv.DictReader(allele_cutoff_file)
            for row in csv_reader:
                cutoffs[row['allele']] = row['allele_specific_cutoff']
        return cutoffs

    @classmethod
    def print_all_allele_cutoffs(cls):
        if not cls.allele_cutoff_dict:
            cls.allele_cutoff_dict = cls.parse_allele_cutoff_file()
        for allele, cutoff in sorted(cls.allele_cutoff_dict.items()):
            print("%s\t%s" % (allele, cutoff))

    @classmethod
    def cutoff_for_allele(cls, allele):
        if not cls.allele_cutoff_dict:
            cls.allele_cutoff_dict = cls.parse_allele_cutoff_file()
        return cls.allele_cutoff_dict.get(allele, None)

    @abstractmethod
    def valid_allele_names(self):
        pass

    @property
    @abstractmethod
    def needs_epitope_length(self):
        pass

    def check_allele_valid(self, allele):
        valid_alleles = self.valid_allele_names()
        if allele not in valid_alleles:
            sys.exit("Allele %s not valid for method %s. Run `pvacseq valid_alleles %s` for a list of valid allele names." % (allele, self.__class__.__name__, self.__class__.__name__))

class MHCI(PredictionClass, metaclass=ABCMeta):
    @property
    def needs_epitope_length(self):
        return True

class MHCflurry(MHCI):
    def valid_allele_names(self):
        predictor = Class1AffinityPredictor.load()
        return predictor.supported_alleles

    def check_length_valid_for_allele(self, length, allele):
        return True

    def valid_lengths_for_allele(self, allele):
        return [8,9,10,11,12,13,14]

    def determine_neoepitopes(self, sequence, length):
        epitopes = []
        for i in range(0, len(sequence)-length+1):
            epitopes.append(sequence[i:i+length])
        return epitopes

    def predict(self, input_file, allele, epitope_length, iedb_executable_path, iedb_retries):
        predictor = Class1AffinityPredictor.load()
        results = pd.DataFrame()
        for line in input_file:
            match = re.search('^>([0-9]+)$', line)
            if match:
                seq_num = match.group(1)
            else:
                epitopes = self.determine_neoepitopes(line.rstrip(), epitope_length)
                if len(epitopes) > 0:
                    df = predictor.predict_to_dataframe(allele=allele, peptides=epitopes)
                    df['seq_num'] = seq_num
                    df['start'] = df.index+1
                    df.rename(columns={'prediction': 'ic50', 'prediction_percentile': 'percentile'}, inplace=True)
                    results = results.append(df)
        return (results, 'pandas')

class MHCnuggetsI(MHCI):
    def valid_allele_names(self):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        alleles_dir       = os.path.join(base_dir, 'tools', 'pvacseq', 'iedb_alleles', 'class_i')
        alleles_file_name = os.path.join(alleles_dir, "MHCnuggets.txt")
        with open(alleles_file_name, 'r') as fh:
            return list(filter(None, fh.read().split('\n')))

    def check_length_valid_for_allele(self, length, allele):
        return True

    def valid_lengths_for_allele(self, allele):
        return [8,9,10,11,12,13,14]

    def find_neoepitopes(self, sequence, length):
        epitopes = defaultdict(list)
        for i in range(0, len(sequence)-length+1):
            epitope = sequence[i:i+length]
            epitopes[epitope].append(i+1)
        return epitopes

    def predict(self, input_file, allele, epitope_length, iedb_executable_path, iedb_retries):
        epitope_seq_nums = defaultdict(list)
        for line in input_file:
            match = re.search('^>([0-9]+)$', line)
            if match:
                seq_num = match.group(1)
            else:
                epitopes = self.find_neoepitopes(line.rstrip(), epitope_length)
                for epitope, starts in epitopes.items():
                    for start in starts:
                        epitope_seq_nums[epitope].append((seq_num, start))
        tmp_file = tempfile.NamedTemporaryFile('w', delete=False)
        for epitope in epitope_seq_nums.keys():
            tmp_file.write("{}\n".format(epitope))
        tmp_file.close()
        tmp_output_file = tempfile.NamedTemporaryFile('r', delete=False)
        predict('I', tmp_file.name, allele.replace('*', ''), output=tmp_output_file.name)
        tmp_output_file.close()
        df = pd.read_csv(tmp_output_file.name)
        processed_df = pd.DataFrame()
        for index, row in df.iterrows():
            seq_nums = epitope_seq_nums[row['peptide']]
            for seq_num, start in seq_nums:
                new_row = row.copy()
                new_row['seq_num'] = seq_num
                new_row['start'] = start
                new_row['allele'] = allele
                processed_df = processed_df.append(new_row)
        processed_df['start'] = pd.to_numeric(processed_df['start'], downcast='integer')
        processed_df = processed_df[['peptide', 'ic50', 'seq_num', 'start', 'allele']]
        return (processed_df, 'pandas')

class IEDBMHCI(MHCI, IEDB, metaclass=ABCMeta):
    @property
    def url(self):
        return 'http://tools-cluster-interface.iedb.org/tools_api/mhci/'

    def parse_iedb_allele_file(self):
        #Ultimately we probably want this method to call out to IEDB but their command is currently broken
        #curl --data "method=ann&species=human" http://tools-api.iedb.org/tools_api/mhci/
        base_dir               = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        iedb_alleles_dir       = os.path.join(base_dir, 'tools', 'pvacseq', 'iedb_alleles', 'class_i')
        iedb_alleles_file_name = os.path.join(iedb_alleles_dir, "%s.tsv" % self.iedb_prediction_method)
        alleles = {}
        with open(iedb_alleles_file_name) as iedb_alleles_file:
            tsv_reader = csv.DictReader(iedb_alleles_file, delimiter='\t')
            for row in tsv_reader:
                allele = row['MHC']
                if allele not in alleles.keys():
                    alleles[allele] = []
                alleles[allele].append(int(row['PeptideLength']))
        return alleles

    def valid_allele_names(self):
        method = self.iedb_prediction_method
        if not self.valid_allele_names_dict:
            self.valid_allele_names_dict = self.parse_iedb_allele_file()
        return self.valid_allele_names_dict.keys()

    def valid_lengths_for_allele(self, allele):
        method = self.iedb_prediction_method
        if not self.valid_allele_names_dict:
            self.valid_allele_names_dict = self.parse_iedb_allele_file()
        return self.valid_allele_names_dict[allele]

    def check_length_valid_for_allele(self, length, allele):
        valid_lengths = self.valid_lengths_for_allele(allele)
        if length not in valid_lengths:
            sys.exit("Length %s not valid for allele %s and method %s." % (length, allele, self.iedb_prediction_method))

    def iedb_executable_params(self, iedb_executable_path, method, allele, input_file, epitope_length):
        return "{} {} {} {} {}".format(iedb_executable_path, method, allele, str(epitope_length), input_file.name)

class NetMHC(IEDBMHCI):
    @property
    def iedb_prediction_method(self):
        return 'ann'

class NetMHCpan(IEDBMHCI):
    @property
    def iedb_prediction_method(self):
        return 'netmhcpan'

class SMMPMBEC(IEDBMHCI):
    @property
    def iedb_prediction_method(self):
        return 'smmpmbec'

class SMM(IEDBMHCI):
    @property
    def iedb_prediction_method(self):
        return 'smm'

class NetMHCcons(IEDBMHCI):
    @property
    def iedb_prediction_method(self):
        return 'netmhccons'

class PickPocket(IEDBMHCI):
    @property
    def iedb_prediction_method(self):
        return 'pickpocket'

class MHCII(PredictionClass, metaclass=ABCMeta):
    @property
    def needs_epitope_length(self):
        return False

class MHCnuggetsII(MHCII):
    def valid_allele_names(self):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        alleles_dir       = os.path.join(base_dir, 'tools', 'pvacseq', 'iedb_alleles', 'class_ii')
        alleles_file_name = os.path.join(alleles_dir, "MHCnuggets.txt")
        with open(alleles_file_name, 'r') as fh:
            return list(filter(None, fh.read().split('\n')))

    def check_length_valid_for_allele(self, length, allele):
        return True

    def valid_lengths_for_allele(self, allele):
        return [15]

    def find_neoepitopes(self, sequence):
        epitopes = defaultdict(list)
        for i in range(0, len(sequence)-14):
            epitope = sequence[i:i+15]
            epitopes[epitope].append(i+1)
        return epitopes

    def predict(self, input_file, allele, epitope_length, iedb_executable_path, iedb_retries):
        epitope_seq_nums = defaultdict(list)
        for line in input_file:
            match = re.search('^>([0-9]+)$', line)
            if match:
                seq_num = match.group(1)
            else:
                epitopes = self.find_neoepitopes(line.rstrip())
                for epitope, starts in epitopes.items():
                    for start in starts:
                        epitope_seq_nums[epitope].append((seq_num, start))
        tmp_file = tempfile.NamedTemporaryFile('w', delete=False)
        for epitope in epitope_seq_nums.keys():
            tmp_file.write("{}\n".format(epitope))
        tmp_file.close()
        tmp_output_file = tempfile.NamedTemporaryFile('r', delete=False)
        mhcnuggets_allele = "HLA-{}".format(allele).replace('*', '')
        predict('II', tmp_file.name, mhcnuggets_allele, output=tmp_output_file.name)
        tmp_output_file.close()
        df = pd.read_csv(tmp_output_file.name)
        processed_df = pd.DataFrame()
        for index, row in df.iterrows():
            seq_nums = epitope_seq_nums[row['peptide']]
            for seq_num, start in seq_nums:
                new_row = row.copy()
                new_row['seq_num'] = seq_num
                new_row['start'] = start
                new_row['allele'] = allele
                processed_df = processed_df.append(new_row)
        processed_df['start'] = pd.to_numeric(processed_df['start'], downcast='integer')
        processed_df = processed_df[['peptide', 'ic50', 'seq_num', 'start', 'allele']]
        return (processed_df, 'pandas')

class IEDBMHCII(MHCII, IEDB, metaclass=ABCMeta):
    @property
    def url(self):
        return 'http://tools-cluster-interface.iedb.org/tools_api/mhcii/'

    def parse_iedb_allele_file(self):
        #Ultimately we probably want this method to call out to IEDB but their command is currently broken
        #curl --data "method=ann&species=human" http://tools-api.iedb.org/tools_api/mhci/
        base_dir               = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        iedb_alleles_dir       = os.path.join(base_dir, 'tools', 'pvacseq', 'iedb_alleles', 'class_ii')
        iedb_alleles_file_name = os.path.join(iedb_alleles_dir, "%s.tsv" % self.iedb_prediction_method)
        alleles = []
        with open(iedb_alleles_file_name) as iedb_alleles_file:
            for row in iedb_alleles_file:
                alleles.append(row.rstrip())
        return alleles

    def valid_allele_names(self):
        method = self.iedb_prediction_method
        if not self.valid_allele_names_dict:
            self.valid_allele_names_dict = self.parse_iedb_allele_file()
        return self.valid_allele_names_dict

    def iedb_executable_params(self, iedb_executable_path, method, allele, input_file, epitope_length):
        allele = allele.replace('-DPB', '/DPB').replace('-DQB', '/DQB')
        return "{} {} {} {}".format(iedb_executable_path, method, allele, input_file.name)

class NetMHCIIpan(IEDBMHCII):
    @property
    def iedb_prediction_method(self):
        return 'NetMHCIIpan'

class NNalign(IEDBMHCII):
    @property
    def iedb_prediction_method(self):
        return 'nn_align'

class SMMalign(IEDBMHCII):
    @property
    def iedb_prediction_method(self):
        return 'smm_align'
