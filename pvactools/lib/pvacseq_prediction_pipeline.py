import os

from pvactools.lib.prediction_pipeline import PredictionPipeline
from pvactools.lib.output_parser import PvacseqOutputParser
from pvactools.tools.pvacseq.generate_protein_fasta import PvacseqGenerateProteinFasta

class PvacseqPredictionPipeline(PredictionPipeline):
    def call_parser(self, parser_arguments):
        PvacseqOutputParser(**parser_arguments).execute()

    def generate_fasta(self):
        fasta_file = os.path.join(self.output_dir, "{}.fasta".format(self.sample_name))
        params = {
            'transcripts_fasta': self.transcript_fasta,
            'output_file': fasta_file,
            'flanking_sequence_length': max(self.epitope_lengths) - 1,
            'mutant_only': False,
        }
        PvacseqGenerateProteinFasta(**params).execute()
