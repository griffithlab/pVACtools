import csv
import argparse
import re
import os
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
import pandas as pd
from abc import ABCMeta, abstractmethod

from pvactools.lib.run_utils import *
import pvactools.lib.sort
from pvactools.lib.prediction_class import PredictionClass

class TopScoreFilter(metaclass=ABCMeta):
    @classmethod
    def parser(cls, tool):
        parser = argparse.ArgumentParser(
            '%s top_score_filter' % tool,
            description="Pick the best neoepitope for each variant",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            'input_file',
            help="The final report .tsv file to filter."
        )
        parser.add_argument(
            'output_file',
            help="Output .tsv file containing only the list of the top "
                 + "epitope per variant."
        )
        parser.add_argument(
            '-m', '--top-score-metric',
            choices=['lowest', 'median'],
            default='median',
            help="The ic50 scoring metric to use for filtering. "
                 + "lowest: Use the best MT Score (i.e. the lowest MT ic50 binding score of all chosen prediction methods). "
                 + "median: Use the median MT Score (i.e. the median MT ic50 binding score of all chosen prediction methods)."
        )
        if tool == 'pvacseq':
            parser.add_argument(
                "--maximum-transcript-support-level", type=int,
                help="When determining the top peptide, only consider those entries that meet this threshold for the Ensembl transcript support level (TSL). "
                     + "Transcript support level needs to be <= this cutoff to be considered.",
                default=1,
                choices=[1,2,3,4,5]
            )
            parser.add_argument(
                '-b', '--binding-threshold', type=int,
                help="When determining the top peptide, only peptides passing the anchor criteria are considered. This criteria is failed if "
                     + "all mutated amino acids of a peptide (Pos) are at an anchor position and the WT peptide has good binding (IC50 WT < binding_threshold).",
                default=500
            )
            parser.add_argument(
                '--allele-specific-binding-thresholds',
                help="When determining the top peptide and evaluating the anchor criteria, use allele-specific binding thresholds. "
                     + "If an allele does not have a special threshold value, the `--binding-threshold` value will be used.",
                default=False,
                action='store_true',
            )

            parser.add_argument(
                "--allele-specific-anchors",
                help="When determining the top peptide and evaluating the anchor criteria, use allele-specific anchor positions. This option "
                     + "is available for 8, 9, 10, and 11mers and only for HLA-A, B, and C alleles. If this option is "
                     + "not enabled or as a fallback for unsupported lengths and alleles, the default positions of 1, "
                     + "2, epitope length - 1, and epitope length are used. Please see https://doi.org/10.1101/2020.12.08.416271 "
                     + "for more details.",
                default=False,
                action='store_true',
            )
            parser.add_argument(
                "--anchor-contribution-threshold", type=float_range(0.5,0.9),
                help="For determining the top peptide and evaluating the anchor criteria using allele-specific anchors, each position is assigned a score based on how binding is "
                     + "influenced by mutations. From these scores, the relative contribution of each position to the "
                     + "overall binding is calculated. Starting with the highest relative contribution, positions whose "
                     + "scores together account for the selected contribution threshold are assigned as anchor locations. "
                     + " As a result, a higher threshold leads to the inclusion of more positions to be considered anchors.",
                default=0.8
            )
        return parser

    @abstractmethod
    def execute(self):
        raise Exception("Must implement method in child class")

    @abstractmethod
    def find_best_line(self, lines):
        raise Exception("Must implement method in child class")


class PvacseqTopScoreFilter(TopScoreFilter, metaclass=ABCMeta):
    def __init__(
        self,
        input_file,
        output_file, 
        top_score_metric="median",
        binding_threshold=500,
        allele_specific_binding_thresholds=False,
        maximum_transcript_support_level=1,
        allele_specific_anchors=False,
        anchor_contribution_threshold=0.8
    ):
        self.input_file = input_file
        self.output_file = output_file
        self.top_score_metric = top_score_metric
        if self.top_score_metric == 'median':
            self.mt_top_score_metric = "Median"
            self.wt_top_score_metric = "Median"
        else:
            self.mt_top_score_metric = "Best"
            self.wt_top_score_metric = "Corresponding"
        self.binding_threshold = binding_threshold
        self.use_allele_specific_binding_thresholds = allele_specific_binding_thresholds
        self.hla_types = pd.read_csv(self.input_file, delimiter="\t", usecols=["HLA Allele"])['HLA Allele'].unique()
        allele_specific_binding_thresholds = {}
        for hla_type in self.hla_types:
            threshold = PredictionClass.cutoff_for_allele(hla_type)
            if threshold is None:
                allele_specific_binding_thresholds[hla_type] = self.binding_threshold
            else:
                allele_specific_binding_thresholds[hla_type] = float(threshold)
        self.allele_specific_binding_thresholds = allele_specific_binding_thresholds
        self.maximum_transcript_support_level = maximum_transcript_support_level
        self.allele_specific_anchors = allele_specific_anchors
        self.anchor_contribution_threshold = anchor_contribution_threshold
        anchor_probabilities = {}
        for length in [8, 9, 10, 11]:
            base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
            file_name = os.path.join(base_dir, 'tools', 'pvacview', 'data', "Normalized_anchor_predictions_{}_mer.tsv".format(length))
            probs = {}
            with open(file_name, 'r') as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for line in reader:
                    hla = line.pop('HLA')
                    probs[hla] = line
            anchor_probabilities[length] = probs
        self.anchor_probabilities = anchor_probabilities

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter = "\t")
            writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames = reader.fieldnames, extrasaction = 'ignore')
            writer.writeheader()
            lines_per_variant = defaultdict(list)
            for line in reader:
                if line["{} MT IC50 Score".format(self.mt_top_score_metric)] == 'NA':
                    continue
                chromosome = line['Chromosome']
                start = line['Start']
                stop = line['Stop']
                ref = line['Reference']
                var = line['Variant']
                index = '%s.%s.%s.%s.%s' % (chromosome, start, stop, ref, var)
                lines_per_variant[index].append(line)

            filtered_lines = []
            for index, lines in lines_per_variant.items():
                lines = sorted(lines, key = itemgetter('Transcript'))
                variant_tracker = defaultdict(set)
                for line in lines:
                    variant, consequence = line['Index'].split('.', 1)
                    variant_tracker[consequence].add(variant)
                transcripts_with_same_epitopes = defaultdict(list)
                for transcript, transcript_lines in groupby(lines, key = itemgetter('Transcript')):
                    transcript_lines = list(transcript_lines)
                    epitopes = ','.join(sorted([x['MT Epitope Seq'] for x in transcript_lines]))
                    transcripts_with_same_epitopes[epitopes].append(transcript)
                for transcripts in transcripts_with_same_epitopes.values():
                    transcript_set_lines = [x for x in lines if x['Transcript'] in transcripts]
                    best_line = self.find_best_line(transcript_set_lines)
                    filtered_lines.append(best_line)
                    best_line_variant, best_line_consequence = best_line['Index'].split('.', 1)
                    for variant in variant_tracker[best_line_consequence]:
                        if variant != best_line_variant:
                            duplicate_variant_line = best_line.copy()
                            duplicate_variant_line['Index'] = "{}.{}".format(variant, best_line_consequence)
                            filtered_lines.append(duplicate_variant_line)


            sorted_rows = pvactools.lib.sort.default_sort(filtered_lines, self.top_score_metric)
            writer.writerows(sorted_rows)

    def find_best_line(self, lines):
        #get all entries with Biotype 'protein_coding'
        biotype_lines = [x for x in lines if x['Biotype'] == 'protein_coding']
        #if there are none, reset to previous dataset
        if len(biotype_lines) == 0:
            biotype_lines = lines

        #subset protein_coding dataset to only include entries with a TSL < maximum_transcript_support_level
        tsl_lines = [x for x in biotype_lines if x['Transcript Support Level'] != 'NA' and x['Transcript Support Level'] != 'Not Supported' and int(x['Transcript Support Level']) < self.maximum_transcript_support_level]
        #if this results in an empty dataset, reset to previous dataset
        if len(tsl_lines) == 0:
            tsl_lines = biotype_lines

        #subset tsl dataset to only include entries with no problematic positions
        if 'Problematic Positions' in lines[0]:
            prob_pos_lines = [x for x in tsl_lines if x['Problematic Positions'] == "None"]
            #if this results in an empty dataset, reset to previous dataset
            if len(prob_pos_lines) == 0:
                prob_pos_lines = tsl_lines
        else:
            prob_pos_lines = tsl_lines

        #subset prob_pos dataset to only include entries that pass the anchor position check
        anchor_residue_pass_lines = [x for x in prob_pos_lines if self.is_anchor_residue_pass(x)]
        if len(anchor_residue_pass_lines) == 0:
            anchor_residue_pass_lines = prob_pos_lines

        #determine the entry with the lowest IC50 Score, lowest TSL, and longest Transcript
        tsl_sort_criteria = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, 'NA': 6, 'Not Supported': 6}
        for line in anchor_residue_pass_lines:
            line['TSL Sort'] = tsl_sort_criteria[line['Transcript Support Level']]
        sorted_anchor_residue_pass_lines = sorted(anchor_residue_pass_lines, key=lambda d: (float(d["{} MT IC50 Score".format(self.mt_top_score_metric)]), d['TSL Sort'], -int(d['Transcript Length'])))
        return sorted_anchor_residue_pass_lines[0]

    def get_anchor_positions(self, hla_allele, epitope_length):
        if self.allele_specific_anchors and epitope_length in self.anchor_probabilities and hla_allele in self.anchor_probabilities[epitope_length]:
            probs = self.anchor_probabilities[epitope_length][hla_allele]
            positions = []
            total_prob = 0
            for (pos, prob) in sorted(probs.items(), key=lambda x: x[1], reverse=True):
                total_prob += float(prob)
                positions.append(pos)
                if total_prob > self.anchor_contribution_threshold:
                    return positions

        return [1, 2, epitope_length - 1 , epitope_length]


    def is_anchor_residue_pass(self, line):
        if self.use_allele_specific_binding_thresholds:
            binding_threshold = self.allele_specific_binding_thresholds[line['HLA Allele']]
        else:
            binding_threshold = self.binding_threshold

        anchor_residue_pass = True
        anchors = self.get_anchor_positions(line['HLA Allele'], len(line['MT Epitope Seq']))
        # parse out mutation position from str
        position = line["Mutation Position"]
        if '-' in position:
            d_ind = position.index('-')
            if all(pos in anchors for pos in range(int(position[0:d_ind]), int(position[d_ind+1:])+1)):
                if line["{} WT IC50 Score".format(self.wt_top_score_metric)] == "NA":
                    anchor_residue_pass = False
                elif float(line["{} WT IC50 Score".format(self.wt_top_score_metric)]) < binding_threshold:
                    anchor_residue_pass = False
        elif position != "NA":
            if int(float(position)) in anchors:
                if line["{} WT IC50 Score".format(self.wt_top_score_metric)] == "NA":
                    anchor_residue_pass = False
                elif float(line["{} WT IC50 Score".format(self.wt_top_score_metric)]) < binding_threshold:
                    anchor_residue_pass = False
        return anchor_residue_pass



class PvacfuseTopScoreFilter(TopScoreFilter, metaclass=ABCMeta):
    def __init__(self, input_file, output_file, top_score_metric="median"):
        self.input_file = input_file
        self.output_file = output_file
        self.top_score_metric = top_score_metric
        if self.top_score_metric == 'median':
            self.formatted_top_score_metric = "Median"
        else:
            self.formatted_top_score_metric = "Best"

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter = "\t")
            writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames = reader.fieldnames, extrasaction = 'ignore')
            writer.writeheader()
            lines_per_variant = defaultdict(list)
            for line in reader:
                if line["{} IC50 Score".format(self.formatted_top_score_metric)] == 'NA':
                    continue
                chromosome = line['Chromosome']
                start = line['Start']
                stop = line['Stop']
                index = '%s.%s.%s' % (chromosome, start, stop)
                lines_per_variant[index].append(line)

            filtered_lines = []
            for index, lines in lines_per_variant.items():
                lines = sorted(lines, key = itemgetter('Mutation'))
                variant_tracker = defaultdict(set)
                for line in lines:
                    variant, consequence = line['Mutation'].split('.', 1)
                    variant_tracker[consequence].add(variant)
                transcripts_with_same_epitopes = defaultdict(list)
                for transcript, transcript_lines in groupby(lines, key = itemgetter('Mutation')):
                    transcript_lines = list(transcript_lines)
                    epitopes = ','.join(sorted([x['Epitope Seq'] for x in transcript_lines]))
                    transcripts_with_same_epitopes[epitopes].append(transcript)
                for transcripts in transcripts_with_same_epitopes.values():
                    transcript_set_lines = [x for x in lines if x['Mutation'] in transcripts]
                    best_line = self.find_best_line(transcript_set_lines)
                    filtered_lines.append(best_line)
                    best_line_variant, best_line_consequence = best_line['Mutation'].split('.', 1)
                    for variant in variant_tracker[best_line_consequence]:
                        if variant != best_line_variant:
                            duplicate_variant_line = best_line.copy()
                            duplicate_variant_line['Mutation'] = "{}.{}".format(variant, best_line_consequence)
                            filtered_lines.append(duplicate_variant_line)

            sorted_rows = pvactools.lib.sort.pvacbind_sort(filtered_lines, self.top_score_metric)
            writer.writerows(sorted_rows)

    def find_best_line(self, lines):
        for line in lines:
            if line['Expression'] == 'NA':
                line['Expression Sort'] = 0
            else:
                line['Expression Sort'] = float(line['Expression'])
        sorted_lines = sorted(lines, key=lambda d: (float(d["{} IC50 Score".format(self.formatted_top_score_metric)]), -d['Expression Sort']))
        return sorted_lines[0]

class PvacbindTopScoreFilter(TopScoreFilter, metaclass=ABCMeta):
    def __init__(self, input_file, output_file, top_score_metric="median"):
        self.input_file = input_file
        self.output_file = output_file
        self.top_score_metric = top_score_metric
        if self.top_score_metric == 'median':
            self.formatted_top_score_metric = "Median"
        else:
            self.formatted_top_score_metric = "Best"

    def execute(self):
        with open(self.input_file) as input_fh, open(self.output_file, 'w') as output_fh:
            reader = csv.DictReader(input_fh, delimiter = "\t")
            writer = csv.DictWriter(output_fh, delimiter = "\t", fieldnames = reader.fieldnames)
            writer.writeheader()
            lines_per_variant = defaultdict(list)
            for line in reader:
                if line["{} IC50 Score".format(self.formatted_top_score_metric)] == 'NA':
                    continue
                lines_per_variant[line['Mutation']].append(line)

            filtered_lines = []
            for index, lines in lines_per_variant.items():
                best_line = self.find_best_line(lines)
                filtered_lines.append(best_line)

            sorted_rows = pvactools.lib.sort.pvacbind_sort(filtered_lines, self.top_score_metric)
            writer.writerows(sorted_rows)

    def find_best_line(self, lines):
        sorted_lines = sorted(lines, key=lambda d: (float(d["{} IC50 Score".format(self.formatted_top_score_metric)])))
        return sorted_lines[0]
