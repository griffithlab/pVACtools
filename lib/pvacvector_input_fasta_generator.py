import os
import csv

class PvacvectorInputFastaGenerator():
    def __init__(self, pvacseq_tsv, input_vcf, output_dir, n_mer):
        self.input_tsv = pvacseq_tsv
        self.input_vcf = input_vcf
        self.output_dir = output_dir
        self.output_file = os.path.join(self.output_dir, "vector_input.fa")
        self.n_mer = int(n_mer)

    def parse_choosen_epitopes(self):
        mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs = [],[],[],[],[],[]
        with open(self.input_tsv, 'r') as input_f:
            reader = csv.DictReader(input_f, delimiter = "\t")
            for line in reader:
                mut_type = line['Variant Type']
                mutation = line['Mutation']
                pos = line['Protein Position']
                gene_name = line['Gene Name']
                mt_epitope_seq = line['MT Epitope Seq']
                wt_epitope_seq = line['WT Epitope Seq']

                mutations.append(mutation)
                mut_types.append(mut_type)
                (old_AA, new_AA) = mutation.split("/")
                #if position presented as a range, use higher end of range
                if "-" in pos:
                    pos = pos.split("-")
                    pos = pos[1]
                    mut_ID = ("MT." + gene_name + "." +  pos + "fs")
                elif mut_type == "FS":
                    mut_ID = "MT." + gene_name + "." + old_AA + pos + "fs"
                elif mut_type == "missense": 
                    mut_ID = "MT." + gene_name + "."  + old_AA + pos + new_AA
                mut_IDs.append(mut_ID)
                mt_epitope_seqs.append(mt_epitope_seq)
                wt_epitope_seqs.append(wt_epitope_seq)
                transcript_IDs.append(line['Transcript'])
        return mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs

    #get necessary data from initial pvacseq input vcf
    def parse_original_vcf(self):
        with open(self.input_vcf, 'r') as input_f:
            transcripts_dict = {}
            for line in input_f:
                attributes = []
                if line[0] != "#":
                    fields = line.split("\t")
                    info = fields[7]
                    info = info.split("|")
                    transcript_ID, downstr_seq, len_change, full_seq = info[6], info[23], info[24], info[25]
                    attributes.append(full_seq)
                    attributes.append(downstr_seq)
                    attributes.append(len_change)
                    transcripts_dict[transcript_ID] = attributes
        input_f.close()
        return(transcripts_dict)

    def edit_full_seq(self, i, mut_types, mutations, wt_epitope_seqs, mt_epitope_seqs, sub_seq, full_seq, transcripts_dict, transcript_IDs):
        if mut_types[i] == "FS":
            downstr_seq, len_change = transcripts_dict[transcript_IDs[i]][1], int(transcripts_dict[transcript_IDs[i]][2])
            parts = mutations[i].split("/")
            initial = parts[0]
            final = parts[1]

            #handle -/X mutations by appending downstr_seq to next position, instead of overwriting last position
            if initial == "-":
                new_end_of_full_seq = len(full_seq) + len_change - len(downstr_seq)
            #overwrites last position of seq with first position of
            #predicted downstr seq
            else:
                new_end_of_full_seq = len(full_seq) + len_change - len(downstr_seq) - 1
            full_seq = full_seq[:new_end_of_full_seq]
        #handles ex: L/LX mutations by adding sequence that is preserved
        #before the downstr predicted sequence
            if len(final) > 1:
                final = final.replace("X", "")
                full_seq = full_seq + final
            full_seq = full_seq + downstr_seq
        elif mut_types[i] == "missense":
            full_seq = full_seq.replace(wt_epitope_seqs[i], mt_epitope_seqs[i])
        else:
            sys.exit("Mutation not yet handled by this parser")
        return(full_seq)

    #get flanking peptides for the epitope chosen
    def get_sub_seq(self, full_seq, mt_seq):
        beginning = full_seq.find(mt_seq)
        if beginning == -1:
            sys.exit("Error: could not find mutant epitope sequence in mutant full sequence")
        length = len(mt_seq)
        end = beginning + length
        #if eptitope sequence is too close to the beginning or end to get the
        #right amount of flanking peptides, get appropriate length from solely
        #ahead or behind
        len_needed = self.n_mer - length
        if len_needed % 2 != 0:
            front = int(beginning - len_needed / 2)
            back = int(end + len_needed / 2)
        else:
            front = int(beginning - len_needed / 2)
            back = int(end + len_needed / 2)
        if front < 0:
            sub_seq = full_seq[beginning:(beginning + self.n_mer)]
        elif back > len(full_seq):
            sub_seq = full_seq[(end - self.n_mer):end]
        else:
            sub_seq = full_seq[front:back]
        return(sub_seq)

    def write_output_fasta(self, mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs, transcripts_dict):
        with open(self.output_file, 'w') as out_f:
            sub_seq = ""
            full_seq = ""
            for i in range(len(transcript_IDs)):
                full_seq = (transcripts_dict[transcript_IDs[i]])[0] 

                full_seq = self.edit_full_seq(i, mut_types, mutations, wt_epitope_seqs, mt_epitope_seqs, sub_seq, full_seq, transcripts_dict, transcript_IDs)

                sub_seq = self.get_sub_seq(full_seq, mt_epitope_seqs[i])
                out_f.write(">" + mut_IDs[i] + "\n")
                out_f.write(sub_seq + "\n")
                print("ID: " + mut_IDs[i] + ", sequence: " + sub_seq)
        out_f.close()
        print("FASTA file written")
        return()


    def execute(self):
        mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs = self.parse_choosen_epitopes()
        transcripts_dict = self.parse_original_vcf()
        self.write_output_fasta(mut_IDs, mutations, mut_types, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs, transcripts_dict)
