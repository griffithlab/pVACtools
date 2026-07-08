from pvactools.lib.prediction_pipeline import PredictionPipeline
from pvactools.lib.output_parser import PvacbindOutputParser

class PvacbindPredictionPipeline(PredictionPipeline):
    def call_parser(self, parser_arguments):
        PvacbindOutputParser(**parser_arguments).execute()

    def generate_fasta(self):
        pass
