from pvactools.tools.pvaccompare.run_utils import *
import logging


class CompareAggregatedTSV:
    def __init__(self, input_file1, input_file2, columns_to_compare):
        self.input_file1 = input_file1
        self.input_file2 = input_file2
        self.df1, self.df2 = load_tsv_files(self.input_file1, self.input_file2)
        self.contains_id = True
        self.replaced_id = False
        self.ID_replacement_cols = ["Gene", "AA Change"]
        self.columns_to_compare = columns_to_compare
        self.original_id_columns = [
            "Chromosome",
            "Start",
            "Stop",
            "Reference",
            "Variant",
        ]

    def check_id(self):
        """
        Purpose:    Replace ID with Gene-AA_change if needed
        Modifies:   self.contains_id, self.replaced_id
        Returns:    None
        """
        if "ID" not in self.df1.columns or "ID" not in self.df2.columns:
            self.contains_id = False
            can_replace = True
            for col in self.ID_replacement_cols:
                if col not in self.df1.columns or col not in self.df2.columns:
                    can_replace = False
            if can_replace:
                self.combine_gene_and_AA_change()
                logging.info("\u2022 Replaced ID with Gene and AA Change")
                self.replaced_id = True

    def combine_gene_and_AA_change(self):
        """
        Purpose:    Combines Gene and AA_Change into a singular unique ID column in both dataframes
        Modifies:   df1 and df2
        Returns:    None
        """
        self.df1["ID"] = (
            self.df1[self.ID_replacement_cols[0]].astype(str)
            + " ("
            + self.df1[self.ID_replacement_cols[1]].astype(str)
            + ")"
        )
        self.df2["ID"] = (
            self.df2[self.ID_replacement_cols[0]].astype(str)
            + " ("
            + self.df2[self.ID_replacement_cols[1]].astype(str)
            + ")"
        )

    def split_id_columns(self):
        """
        Purpose:    Split the combined ID column into Chromosome, Start, Stop, Reference, and Variant
        Modifies:   df1 and df2
        Returns:    None
        """

        def split_id(df):
            split_df = df["ID"].str.split("-", expand=True)
            split_df.columns = self.original_id_columns
            for col in self.original_id_columns:
                df[col] = split_df[col]

        split_id(self.df1)
        split_id(self.df2)
