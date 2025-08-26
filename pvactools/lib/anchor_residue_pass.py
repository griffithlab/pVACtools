import pandas as pd
import os
import csv
import ast

class AnchorResiduePass:
    def __init__(self, binding_threshold, use_allele_specific_binding_thresholds, allele_specific_binding_thresholds, allele_specific_anchors, anchor_contribution_threshold, wt_top_score_metric=None):
        self.default_binding_threshold = binding_threshold
        self.use_allele_specific_binding_thresholds = use_allele_specific_binding_thresholds
        self.allele_specific_binding_thresholds = allele_specific_binding_thresholds
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

        mouse_anchor_positions = {}
        for length in [8, 9, 10, 11]:
            base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
            file_name = os.path.join(base_dir, 'tools', 'pvacview', 'data', "mouse_anchor_predictions_{}_mer.tsv".format(length))
            values = {}
            with open(file_name, 'r') as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for line in reader:
                    allele = line.pop('Allele')
                    values[allele] = {int(k): ast.literal_eval(v) for k, v in line.items()}
            mouse_anchor_positions[length] = values
        self.mouse_anchor_positions = mouse_anchor_positions

        self.use_allele_specific_anchors = allele_specific_anchors
        self.anchor_contribution_threshold = anchor_contribution_threshold
        self.wt_top_score_metric = wt_top_score_metric

    def binding_threshold(self, allele):
        if self.use_allele_specific_binding_thresholds and allele in self.allele_specific_binding_thresholds:
            binding_threshold = self.allele_specific_binding_thresholds[allele]
        else:
            binding_threshold = self.default_binding_threshold
        return binding_threshold

    def get_anchor_positions(self, hla_allele, epitope_length):
        if self.use_allele_specific_anchors and epitope_length in self.anchor_probabilities and hla_allele in self.anchor_probabilities[epitope_length]:
            probs = self.anchor_probabilities[epitope_length][hla_allele]
            positions = []
            total_prob = 0
            for (pos, prob) in sorted(probs.items(), key=lambda x: x[1], reverse=True):
                total_prob += float(prob)
                positions.append(int(pos))
                if total_prob > self.anchor_contribution_threshold:
                    return positions
        elif self.use_allele_specific_anchors and epitope_length in self.mouse_anchor_positions and hla_allele in self.mouse_anchor_positions[epitope_length]:
            values = self.mouse_anchor_positions[epitope_length][hla_allele]
            positions = [pos for pos, val in values.items() if val]
            return positions
        return [1, 2, epitope_length - 1 , epitope_length]

    def is_anchor_residue_pass(self, mutation):
        if 'Allele' in mutation:
            allele = mutation['Allele']
            peptide = mutation['Best Peptide']
            position = mutation['Pos']
            wt_score = mutation['IC50 WT']
        else:
            allele = mutation['HLA Allele']
            peptide = mutation['MT Epitope Seq']
            position = mutation['Mutation Position']
            wt_score = mutation["{} WT IC50 Score".format(self.wt_top_score_metric)]
        anchors = self.get_anchor_positions(allele, len(peptide))
        # parse out mutation positions from str
        if position == 'NA' or pd.isna(position):
            return True
        else:
            positions = position.split(", ")
            if len(positions) > 2:
                return True
            anchor_residue_pass = True
            if all(int(pos) in anchors for pos in positions):
                if wt_score == 'NA':
                    anchor_residue_pass = False
                elif float(wt_score) < self.binding_threshold(allele):
                    anchor_residue_pass = False
            return anchor_residue_pass
