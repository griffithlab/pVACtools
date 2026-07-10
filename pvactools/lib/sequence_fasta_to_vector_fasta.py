import os
import json
import logging

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class SequenceFastaToVectorFasta():
    def __init__(self, **kwargs):
        self.input_file         = kwargs['input_file']
        self.output_file        = kwargs['output_file']
        self.epitope_length     = kwargs['epitope_length']
        self.spacer             = kwargs['spacer']
        self.junctions_to_test  = kwargs['junctions_to_test']
        self.clip_length        = kwargs['clip_length']

    def execute(self):
        seq_dict = dict()
        best_peptide = dict()
        for record in SeqIO.parse(self.input_file, "fasta"):
            seq_dict[record.id] = str(record.seq)
            description = record.description.replace("{} ".format(record.id), "")
            if description != "":
               try:
                   best_peptide[record.id] = json.loads(description)['Best Peptide']
               except:
                   pass

        records = []
        wingspan_length = self.epitope_length - 1
        warnings = set()
        for (seq1, seq2) in self.junctions_to_test:
            seq1_seq = seq_dict[seq1]
            seq2_seq = seq_dict[seq2]
            for left_clip_length in range(0, self.clip_length+1):
                for right_clip_length in range(0, self.clip_length+1):
                    #These combinations would've already been tested in previous attempts with lower clip lengths and can be skipped
                    if left_clip_length < self.clip_length and right_clip_length < self.clip_length:
                        continue
                    if seq1 in best_peptide:
                        seq1_best_peptide = best_peptide[seq1]
                        last_position = seq1_seq.rindex(seq1_best_peptide) + len(seq1_best_peptide)
                        end_distance = len(seq1_seq) - last_position
                        if left_clip_length > end_distance:
                            warnings.add(f"Clipping {left_clip_length} amino acids off the end of peptide {seq1} would clip the best peptide. Skipping.")
                            continue
                    if seq2 in best_peptide:
                        seq2_best_peptide = best_peptide[seq2]
                        first_position = seq2_seq.index(seq2_best_peptide)
                        if right_clip_length > first_position:
                            warnings.add(f"Clipping {right_clip_length} amino acids off the start of peptide {seq2} would clip the best peptide. Skipping.")
                            continue
                    trunc_seq1 = seq1_seq[(len(seq1_seq) - wingspan_length - left_clip_length):(len(seq1_seq) - left_clip_length)]
                    trunc_seq2 = seq2_seq[(0 + right_clip_length):wingspan_length + right_clip_length]

                    if self.spacer != 'None':
                        seq_ID = f"{seq1}|{left_clip_length}|{self.spacer}|{right_clip_length}|{seq2}"
                        records.append(SeqRecord(Seq(trunc_seq1 + self.spacer + trunc_seq2), id=seq_ID, description=""))
                    else:
                        seq_ID = f"{seq1}|{left_clip_length}|{right_clip_length}|{seq2}"
                        records.append(SeqRecord(Seq(trunc_seq1 + trunc_seq2), id=seq_ID, description=""))
        for warning in list(warnings):
            logging.info(warning)

        SeqIO.write(records, self.output_file, "fasta")
