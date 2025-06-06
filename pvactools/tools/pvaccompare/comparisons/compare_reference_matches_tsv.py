from pvactools.tools.pvaccompare.run_utils import *
import logging


class CompareReferenceMatchesTSV:
    def __init__(self, input_file1, input_file2, columns_to_compare):
        self.input_file1 = input_file1
        self.input_file2 = input_file2
        self.df1, self.df2 = load_tsv_files(self.input_file1, self.input_file2)
        self.columns_to_compare = columns_to_compare
        self.run_notes = []
        self.hits_file1 = {}
        self.hits_file2 = {}
        self.original_id_columns = [
            "Chromosome",
            "Start",
            "Stop",
            "Reference",
            "Variant",
            "Transcript",
            "MT Epitope Seq",
            "Hit ID",
            "Match Start",
            "Match Stop",
        ]

    def check_duplicate_ids(self):
        """
        Purpose:    Checks if duplicate IDs exist in either dataframe
        Modifies:   Nothing
        Returns:    Boolean value
        """
        self.hits_file1 = self.df1["ID"].value_counts().to_dict()
        self.hits_file2 = self.df2["ID"].value_counts().to_dict()

        max_hits_file1 = max(self.hits_file1.values(), default=0)
        max_hits_file2 = max(self.hits_file2.values(), default=0)

        if max_hits_file1 > 1 or max_hits_file2 > 1:
            if max_hits_file1 > 1 and max_hits_file2 > 1:
                logging.error(
                    "ERROR: Duplicate unique records were found in both files. Writing number of hits only."
                )
                self.run_notes.append(
                    "ERROR: Duplicate unique records were found in both files. Writing number of hits only."
                )
            elif max_hits_file1 > 1:
                logging.error(
                    "ERROR: Duplicate unique records were found in file 1. Writing number of hits only."
                )
                self.run_notes.append(
                    "ERROR: Duplicate unique records were found in file 1. Writing number of hits only."
                )
            else:
                logging.error(
                    "ERROR: Duplicate unique records were found in file 2. Writing number of hits only."
                )
                self.run_notes.append(
                    "ERROR: Duplicate unique records were found in file 2. Writing number of hits only."
                )
            self.hits_file1 = dict(
                sorted(self.hits_file1.items(), key=lambda x: extract_id_parts(x[0]))
            )
            self.hits_file2 = dict(
                sorted(self.hits_file2.items(), key=lambda x: extract_id_parts(x[0]))
            )
            return True
        else:
            return False
