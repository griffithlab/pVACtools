from abc import ABCMeta, abstractmethod
import os
import csv
import sys

class PredictionClass(metaclass=ABCMeta):
    valid_allele_names_dict = {}

    @classmethod
    def prediction_classes(cls):
        prediction_classes = []
        for subclass in cls.__subclasses__():
            for prediction_class in subclass.__subclasses__():
                prediction_classes.append(prediction_class)
        return prediction_classes

    @classmethod
    def prediction_methods(cls):
        return sorted([prediction_class.__name__ for prediction_class in cls.prediction_classes()])

    @classmethod
    def iedb_prediction_methods(cls):
        return [prediction_class().iedb_prediction_method for prediction_class in cls.prediction_classes()]

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
    def prediction_class_for_iedb_prediction_method(cls, method):
        prediction_classes = cls.prediction_classes()
        for prediction_class in prediction_classes:
            prediction_class_object = prediction_class()
            if prediction_class_object.iedb_prediction_method == method:
                return prediction_class_object

    @classmethod
    def prediction_class_name_for_iedb_prediction_method(cls, method):
        return cls.prediction_class_for_iedb_prediction_method(method).__class__.__name__

    @abstractmethod
    def parse_iedb_allele_file(self):
        pass

    @abstractmethod
    def valid_allele_names(self):
        pass

    @property
    @abstractmethod
    def iedb_prediction_method(self):
        pass

    @property
    @abstractmethod
    def url(self):
        pass

    def check_allele_valid(self, allele):
        valid_alleles = self.valid_allele_names()
        if allele not in valid_alleles:
            sys.exit("Allele %s not valid for method %s. Run `pvacseq valid_alleles %s` for a list of valid allele names." % (allele, self.iedb_prediction_method, self.__class__.__name__))

class MHCI(PredictionClass, metaclass=ABCMeta):
    @property
    def url(self):
        return 'http://tools-api.iedb.org/tools_api/mhci/'

    def parse_iedb_allele_file(self):
        #Ultimately we probably want this method to call out to IEDB but their command is currently broken
        #curl --data "method=ann&species=human" http://tools-api.iedb.org/tools_api/mhci/
        base_dir               = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        iedb_alleles_dir       = os.path.join(base_dir, 'iedb_alleles', 'class_i')
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

class NetMHC(MHCI):
    @property
    def iedb_prediction_method(self):
        return 'ann'

class NetMHCpan(MHCI):
    @property
    def iedb_prediction_method(self):
        return 'netmhcpan'

class SMMPMBEC(MHCI):
    @property
    def iedb_prediction_method(self):
        return 'smmpmbec'

class SMM(MHCI):
    @property
    def iedb_prediction_method(self):
        return 'smm'

class NetMHCcons(MHCI):
    @property
    def iedb_prediction_method(self):
        return 'netmhccons'

class PickPocket(MHCI):
    @property
    def iedb_prediction_method(self):
        return 'pickpocket'

class MHCII(PredictionClass, metaclass=ABCMeta):
    @property
    def url(self):
        return 'http://tools-api.iedb.org/tools_api/mhcii/'

    def parse_iedb_allele_file(self):
        #Ultimately we probably want this method to call out to IEDB but their command is currently broken
        #curl --data "method=ann&species=human" http://tools-api.iedb.org/tools_api/mhci/
        base_dir               = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        iedb_alleles_dir       = os.path.join(base_dir, 'iedb_alleles', 'class_ii')
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

class NetMHCIIpan(MHCII):
    @property
    def iedb_prediction_method(self):
        return 'NetMHCIIpan'

class NNalign(MHCII):
    @property
    def iedb_prediction_method(self):
        return 'nn_align'

class SMMalign(MHCII):
    @property
    def iedb_prediction_method(self):
        return 'smm_align'
