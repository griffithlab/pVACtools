from pvactools.tools.pvaccompare.run_utils import *


class CompareUnaggregatedTSV:
    def __init__(self, input_file1, input_file2, columns_to_compare):
        self.input_file1 = input_file1
        self.input_file2 = input_file2
        self.df1, self.df2 = load_tsv_files(self.input_file1, self.input_file2)
        self.columns_to_compare = columns_to_compare

    def create_id_column(self):
        """
        Purpose:    Combines multiple columns into a singular unique ID column in both dataframes
        Modifies:   df1 and df2
        Returns:    None
        """
        id_columns = [
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
        self.df1["ID"] = self.df1[id_columns].apply(
            lambda x: "-".join(map(str, x)), axis=1
        )
        self.df2["ID"] = self.df2[id_columns].apply(
            lambda x: "-".join(map(str, x)), axis=1
        )

        self.df1.drop(columns=id_columns, inplace=True)
        self.df2.drop(columns=id_columns, inplace=True)
