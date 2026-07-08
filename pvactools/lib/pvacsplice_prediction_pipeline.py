import os

from pvactools.lib.prediction_pipeline import PredictionPipeline
from pvactools.lib.output_parser import PvacspliceOutputParser
from pvactools.tools.pvacsplice.generate_protein_fasta import PvacspliceGenerateProteinFasta

class PvacsplicePredictionPipeline(PredictionPipeline):
    def call_parser(self, parser_arguments):
        PvacspliceOutputParser(**parser_arguments).execute()

    def generate_fasta(self):
        fasta_file = os.path.join(self.output_dir, "{}.fasta".format(self.sample_name))
        params = {
            'fasta_file_path': self.transcript_fasta,
            'trimmed_fasta_file_path': fasta_file,
            'flanking_sequence_length': max(self.epitope_lengths) - 1,
            'mutant_only': False,
        }
        PvacspliceGenerateProteinFasta(**params).trim_sequences()
