import os

class BasePipeline:
    @staticmethod
    def file_exists(file_path: str, file_type: str):
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            print(f"{file_type} file already exists. Skipping.")
            exists = True
        else:
            exists = False
        return exists

    def execute(self):
        self.generate_fasta()
        self.fasta_to_kmers()

    def generate_fasta(self):
        self.input_to_tsv()
        self.tsv_to_fasta()

    def choose_final_lengths(self):
        if not self.class_i_hla:
            lengths = self.class_ii_epitope_length
        elif not self.class_ii_hla:
            lengths = self.class_i_epitope_length
        else:
            lengths = self.class_i_epitope_length + self.class_ii_epitope_length
        return lengths
