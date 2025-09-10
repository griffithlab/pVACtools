from pvactools.tools.pvaccompare.run_utils import *


class CompareUnaggregatedTSV:
    def __init__(self, input_file1, input_file2, columns_to_compare):
        self.input_file1 = input_file1
        self.input_file2 = input_file2
        self.df1, self.df2 = load_tsv_files(self.input_file1, self.input_file2)
        self.columns_to_compare = columns_to_compare
        self.original_id_columns = [
            "Chromosome",
            "Start",
            "Stop",
            "Reference",
            "Variant",
            "HLA Allele",
            "Sub-peptide Position",
            "MT Epitope Seq",
            "Index",
        ]
