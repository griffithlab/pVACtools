import os
import re
import pandas as pd
from pyfaidx import Fasta
from Bio.Seq import Seq
from load_ensembl_data import *

class JunctionToFasta():
    def __init__(self, **kwargs):
        self.ref_fasta      = Fasta(kwargs['fasta_path'])
        self.tscript_id     = kwargs['tscript_id']
        self.chrom          = kwargs['chrom']
        self.junction_name  = kwargs['junction_name']
        self.junction_coors = kwargs['junction_coors']
        self.junction_index  = kwargs['index']
        self.anchor         = kwargs['anchor']
        self.strand         = kwargs['strand']
        self.gene_name      = kwargs['gene_name']
        self.output_file    = kwargs['output_file']
        self.output_dir     = kwargs['output_dir']
        if (self.anchor == 'A' and self.strand == 1) or (self.anchor == 'D' and self.strand == -1) or (self.anchor == 'NDA'):
            self.wt_coor  = int(self.junction_coors[1])
            self.alt_coor = int(self.junction_coors[0])
            self.wt_row   = "Genomic coding start"
            self.alt_row  = "Genomic coding end"
            self.reverse  = True
        elif (self.anchor == 'D' and self.strand == 1) or (self.anchor == 'A' and self.strand == -1):
            self.wt_coor  = int(self.junction_coors[0])
            self.alt_coor = int(self.junction_coors[1])
            self.wt_row   = "Genomic coding end"
            self.alt_row  = "Genomic coding start"
            self.reverse  = False
        self.wt_df  = pd.DataFrame()
        self.alt_df = pd.DataFrame()       
        # testing
        # self.test_dir = f'/Users/mrichters/Documents/GitHub/pVACsplice/tests/test_data/junction_to_fasta/{self.tscript_id}'
        # if not os.path.isdir(self.test_dir):
        #     os.mkdir(self.test_dir)

    def create_wt_df(self):
        # load in wt transcript df
        self.wt_df = load_ensembl_data(self.tscript_id)
        # if anchor is D or A, make sure the wt coordinate is inside in the coding region of transcript
        # (alt coordinate won't be present because ensembl only lists the wt coordinates) 
        if self.anchor in ['D', 'A']:
            # look for ref coor index
            index = [index for index,value in self.wt_df[self.wt_row].items() if value == self.wt_coor]
            # value not found in df (Exception)
            if not index:
                print(f'{self.tscript_id} anchor coordinate is not within coding transcript...skipping')
                self.wt_df = pd.DataFrame()
        # if exon skip, check that both coordinates are inside the coding region of transcript
        elif self.anchor == 'NDA':
            # check for presence of both coordinates 
            index_wt = [index for index,value in self.wt_df[self.wt_row].items() if value == self.wt_coor]
            index_alt = [index for index,value in self.wt_df[self.alt_row].items() if value == self.alt_coor]
            if not index_wt or not index_alt:
                # here i can add option to look for next start OR stop codon
                print(f'{self.tscript_id} exon skip: at least 1 junction coordinate is not within coding transcript...skipping')
                self.wt_df = pd.DataFrame()
        # testing
        # self.wt_df.to_csv(f'{self.test_dir}/wt_dataframe.tsv', sep='\t', index=False)    
        return self.wt_df
        

    def create_alt_df(self):
        # copy wt_df so its not directly affected
        self.alt_df = self.wt_df.copy()
        # reverse direction
        if self.reverse == True and (self.anchor != 'NDA'):
            # get rows below the wt coordinate to loop over and find the right index for alt coor insertion
            revised_df = self.alt_df[self.alt_df[self.wt_row] <= self.wt_coor]
            # loop over indexes in revised_df
            # ex: [3,2,1,0] (going to beginning of df)
            index_list = list(reversed(revised_df.index.tolist()))
            for i in index_list:
                # # if i IS in between two wt coordinates and i IS NOT the last index in list (i-1 index exists in df)
                if i-1 in index_list and self.alt_coor < revised_df.loc[i, self.wt_row] and self.alt_coor > revised_df.loc[i-1, self.wt_row]:
                    # create alt_index
                    alt_index = i-1
                    # modify alt_df to include alt coor and return new alt_df
                    self.alt_df.at[alt_index, self.alt_row] = self.alt_coor
                    break
                # if i IS NOT between two wt coordinates and it IS NOT the last index in list
                # go to next index in loop
                elif i != index_list[-1]:
                    continue
                # if i IS NOT between two wt coordinates AND it IS the last index in list
                # return an empty df that will cause an exception in run.py                
                else:
                    # here i can add option to look for next start codon (start lost)
                    print(f'{self.tscript_id} alternate junction coordinate is not within coding transcript...skipping')
                    self.alt_df = pd.DataFrame()
                
        # forward direction
        if self.reverse == False:
            # select rows that are possible places to insert altant coordinate
            revised_df = self.alt_df[self.alt_df[self.wt_row] >= self.wt_coor] 
            # loop over indexes in revised_df
            # ex: [13,14,15] (going to end of df)
            index_list = revised_df.index.tolist()
            for i in index_list:
                # if i IS in between two wt coordinates and i IS NOT the last index in list (i+1 index exists in df)
                if i+1 in index_list and self.alt_coor > revised_df.loc[i, self.wt_row] and self.alt_coor < revised_df.loc[i+1, self.wt_row]:
                    # create alt_index
                    alt_index = i+1
                    # modify alt_df to include alt coor and return new alt_df
                    self.alt_df.at[alt_index, self.alt_row] = self.alt_coor
                    break
                # if i IS NOT between two wt coordinates and it IS NOT the last index in list
                # go to next index in loop
                elif i != index_list[-1]:
                    continue
                # if i IS NOT between two wt coordinates AND it IS the last index in list
                # return an empty df that will cause an exception in run.py
                else:
                    # here i can add option to look for next stop codon (stop lost)
                    print(f'{self.tscript_id} alternate junction coordinate is not within coding transcript...skipping')
                    self.alt_df = pd.DataFrame()    
                
        # exon skip
        # starting knowing both coors are in coding region
        if self.anchor == 'NDA':
            start_index = self.alt_df[self.alt_df[self.wt_row] == self.wt_coor].index.item()
            stop_index = self.alt_df[self.alt_df[self.alt_row] == self.alt_coor].index.item()
            # delete any rows in between
            # [1:] to remove 1st index from index_list (because of range function)
            index_list = list(range(stop_index, start_index))[1:]
            self.alt_df = self.alt_df.drop(index_list)
            
        # testing
        #self.alt_df.to_csv(f'{self.test_dir}/alt_dataframe.tsv', sep='\t', index=False)
        return self.alt_df


    def get_aa_sequence(self, dataframe):
        # pyfaidx has 0-based indexing so subtract 1 from coding exon start positions
        dataframe["Genomic coding start"] = dataframe["Genomic coding start"] -1 
        # create coding_coors column for fasta indexing
        dataframe["coding_coors"] = dataframe["Genomic coding start"].astype(str) + ":" + dataframe["Genomic coding end"].astype(str)
        coordinates = dataframe["coding_coors"].tolist()
        # generate AA sequence from coding exon coordinates (pyfaidx)
        final_seq = []
        for x in coordinates:
            coor = x.replace("'", "")
            seq = f"{self.ref_fasta}['{self.chrom}'][{coor}].seq"
            final_seq.append(seq)
        # string
        str_seq = eval((" + ").join(final_seq))
        # using Seq from Bio.Seq to translate str_seq
        # positive strand
        dna_seq = Seq(str_seq) # creating Seq object
        # negative strand
        if self.strand == -1:
            dna_seq = dna_seq.reverse_complement()
        # make biopython happy by making dna_seq a multiple of 3 
        # adding Ns to end of sequence if remainder != 0
        remainder = len(dna_seq) % 3
        if remainder == 1:
            dna_seq += "NN"
        elif remainder == 2:
            dna_seq += "N"
        aa_seq = str(dna_seq.translate(to_stop=True))
        if aa_seq[0] != 'M':
            print(f'{self.tscript_id} does not begin with start codon...skipping')
            aa_seq = ''

        return aa_seq

    def create_sequence_fasta(self, wt_protein, alt_protein):
        #testing
        #file = f'{self.test_dir}/protein.fasta'
        #os.mkdir(f'{self.working_dir}/results')
        write_str = f'>WT.{self.junction_index}\n{wt_protein}\n>ALT.{self.junction_index}\n{alt_protein}\n'
        if os.path.exists(self.output_file):
            dup_content = re.search(write_str, open(self.output_file, "r").read())
            if dup_content == None:
                with open(self.output_file, "a") as f:
                    f.write(write_str)
        else:
            with open(self.output_file, "a") as f:
                f.write(write_str)



# GBM examples #
# ex 1: KLHL5 A (+, 1) FS # 
# 'ENST00000261425', 'chr4', [39081256, 39081963], 'A', 1, 'KLHL5'

# ex 2: DSG3 D (+, 1) INDEL (1 aa deletion) #
# 'ENST00000257189', 'chr18', [31472788, 31474124], 'D', 1, 'DSG3'

# ex 3: RAB18 NDA (exon skip) (+, 1) FS #
# 'ENST00000356940', 'chr10', [27509930,27532507], 'NDA', 1, 'RAB18'

# ex 4: ATP5G3 A (-, -1) FS #
# 'ENST00000284727', 'chr2', [175178402, 175179136], 'A', -1, 'ATP5G3'

# ex 5: LY6D D (-, -1) FS #
# 'ENST00000301263', 'chr8', [142785451, 142785589], 'D', -1, 'LY6D'

# ex 6: SHPRH NDA (-, -1) stop codon created immediately so no novel peptides possible #
# start 9; stop 11 -remove line 10
# 'ENST00000275233', 'chr6', [145921392, 145922663], 'NDA', -1, 'SHPRH'
