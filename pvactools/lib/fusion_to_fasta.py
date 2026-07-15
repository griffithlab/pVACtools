from abc import ABCMeta
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, translate
from collections import OrderedDict, defaultdict
import csv
import tempfile
import gzip
import shutil

import pvactools.lib.run_utils

class FusionToFasta(metaclass=ABCMeta):
    def __init__(self, **kwargs):
        self.input_file = kwargs['input_file']
        if pvactools.lib.run_utils.is_gz_file(kwargs['transcript_fasta']):
            unzipped_file = tempfile.NamedTemporaryFile('wb')
            with gzip.open(kwargs['transcript_fasta'], "rb") as f_in:
                shutil.copyfileobj(f_in, unzipped_file)
            self.transcript_fasta = unzipped_file.name
        else:
            self.transcript_fasta = kwargs['transcript_fasta']
        self.transcript_fasta_dict_versioned = SeqIO.to_dict(SeqIO.parse(self.transcript_fasta, "fasta"))
        self.transcript_fasta_dict_unversioned = { k.split('.')[0]: v for k, v in self.transcript_fasta_dict_versioned.items() }
        self.downstream_sequence_length = kwargs['downstream_sequence_length']
        self.output_file = kwargs['output_file']

    def execute(self):
        records = []
        with open(self.input_file, 'r') as input_fh:
            reader = csv.DictReader(input_fh, delimiter='\t')
            for line in reader:
                fusion_seq = line['fusion_amino_acid_sequence']
                fusion_pos = int(line['protein_position'])
                fusion_records = []
                five_transcript, three_transcript = line['transcript_name'].split('-')

                full_five_transcript_seq = self.get_transcript_peptide_sequence(five_transcript)
                if full_five_transcript_seq is None:
                    continue
                five_transcript_seq = self.trim_five_transcript_seq(full_five_transcript_seq, fusion_seq, fusion_pos, line['index'])
                if five_transcript_seq is None:
                    continue
                if fusion_seq in five_transcript_seq:
                    continue
                five_transcript_id = f'WT5.{line["index"]}'
                fusion_records.append(SeqRecord(Seq(five_transcript_seq), id=five_transcript_id, description=""))

                if line['variant_type'] == 'inframe_fusion':
                    full_three_transcript_seq = self.get_transcript_peptide_sequence(three_transcript)
                    if full_three_transcript_seq is None:
                        continue
                    three_transcript_seq = self.trim_three_transcript_seq(full_three_transcript_seq, fusion_seq, fusion_pos, line['index'])
                    if three_transcript_seq is None:
                        continue
                    three_transcript_id = f'WT3.{line["index"]}'
                    fusion_records.append(SeqRecord(Seq(three_transcript_seq), id=three_transcript_id, description=""))
                else:
                    if self.downstream_sequence_length is not None:
                        fusion_seq = fusion_seq[:(fusion_pos + self.downstream_sequence_length)]

                fusion_transcript_id = f'MT.{line["index"]}'
                fusion_records.append(SeqRecord(Seq(fusion_seq), id=fusion_transcript_id, description=""))

                records.extend(fusion_records)

        SeqIO.write(records, self.output_file, "fasta")


    def get_transcript_peptide_sequence(self, transcript):
        if '.' in transcript:
            if transcript not in self.transcript_fasta_dict_versioned:
                print(f"{transcript} not found in CDS FASTA file. Skipping")
                return None
            transcript_seq = self.transcript_fasta_dict_versioned[transcript]
        else:
            if transcript not in self.transcript_fasta_dict_unversioned:
                print(f"{transcript} not found in CDS FASTA file. Skipping")
                return None
            transcript_seq = self.transcript_fasta_dict_unversioned[transcript]
        seq = str(translate(transcript_seq.seq))
        if seq.endswith('*'):
            seq = seq[:-1]
        return seq

    def trim_five_transcript_seq(self, full_five_transcript_seq, fusion_seq, fusion_pos, index):
        five_fusion_seq = fusion_seq[:fusion_pos]
        if full_five_transcript_seq.startswith(five_fusion_seq):
            return full_five_transcript_seq
        else:
            if five_fusion_seq in full_five_transcript_seq:
                trimmed_five_transcript_start = full_five_transcript_seq.index(five_fusion_seq)
                return full_five_transcript_seq[trimmed_five_transcript_start:]
            else:
                print(f"{index}: 5' fusion sequence {five_fusion_seq} not found in 5' transcript {full_five_transcript_seq}. Skipping")
                return None

    def trim_three_transcript_seq(self, full_three_transcript_seq, fusion_seq, fusion_pos, index):
        three_fusion_seq = fusion_seq[fusion_pos+1:]
        if full_three_transcript_seq.endswith(three_fusion_seq):
            return full_three_transcript_seq
        else:
            if three_fusion_seq in full_three_transcript_seq:
                trimmed_three_transcript_end = full_three_transcript_seq.index(three_fusion_seq) + len(three_fusion_seq)
                return full_three_transcript_seq[:trimmed_three_transcript_end]
            else:
                print(f"{index}: 3' fusion sequence {three_fusion_seq} not found in 3' transcript {full_three_transcript_seq}. Skipping")
                return None
