import argparse
import sys
import os
import re
import pandas as pd
from pyfaidx import Fasta
from Bio.Seq import Seq
from load_ensembl_data import load_ensembl_data


class JunctionToFasta():
    def __init__(self, fasta_path:str, tscript_id:str, chrom:int, junction_name:str, junction_coors:list, anchor:str, strand:int, gene_name:str, output_dir:str):
        self.wt_df = pd.DataFrame()
        self.mt_df = pd.DataFrame()
        if (anchor == 'A' and strand == 1) or (anchor == 'D' and strand == -1) or (anchor == 'NDA'):
            self.wt_coor = int(junction_coors[1])
            self.mt_coor = int(junction_coors[0])
            self.wt_row = "Genomic coding start"
            self.mt_row = "Genomic coding end"
            self.reverse = True
        elif (anchor == 'D' and strand == 1) or (anchor == 'A' and strand == -1):
            self.wt_coor = int(junction_coors[0])
            self.mt_coor = int(junction_coors[1])
            self.wt_row = "Genomic coding end"
            self.mt_row = "Genomic coding start"
            self.reverse = False
        self.tscript_id = tscript_id
        self.chrom = chrom
        self.junction_name = junction_name
        self.anchor = anchor
        self.strand = strand
        self.gene_name = gene_name
        self.ref_fasta = Fasta(fasta_path)
        self.output_dir = output_dir
        # testing
        # self.test_dir = f'/Users/mrichters/Documents/GitHub/pVACsplice/tests/test_data/junction_to_fasta/{self.tscript_id}'
        # if not os.path.isdir(self.test_dir):
        #     os.mkdir(self.test_dir)

    def createWtDataframe(self):
        # load in wt transcript df
        self.wt_df = load_ensembl_data(self.tscript_id)
        # if anchor is D or A, make sure the wt coordinate is inside in the coding region of transcript
        # (mt coordinate won't be present because ensembl only lists the wt coordinates) 
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
            index_mt = [index for index,value in self.wt_df[self.mt_row].items() if value == self.mt_coor]
            if not index_wt or not index_mt:
                # here i can add option to look for next start OR stop codon
                print(f'{self.tscript_id} exon skip: at least 1 junction coordinate is not within coding transcript...skipping')
                self.wt_df = pd.DataFrame()
        
        # testing
        # self.wt_df.to_csv(f'{self.test_dir}/wt_dataframe.tsv', sep='\t', index=False)    
        
        return self.wt_df
        

    def createMtDataframe(self):
        # copy wt_df so its not directly affected
        self.mut_df = self.wt_df.copy()
        
        # reverse direction
        if self.reverse == True and (self.anchor != 'NDA'):
            # get rows below the wt coordinate to loop over and find the right index for mt coor insertion
            revised_df = self.mut_df[self.mut_df[self.wt_row] <= self.wt_coor]
            # loop over indexes in revised_df
            # ex: [3,2,1,0] (going to beginning of df)
            index_list = list(reversed(revised_df.index.tolist()))
            for i in index_list:
                # # if i IS in between two wt coordinates and i IS NOT the last index in list (i-1 index exists in df)
                if i-1 in index_list and self.mt_coor < revised_df.loc[i, self.wt_row] and self.mt_coor > revised_df.loc[i-1, self.wt_row]:
                    # create mt_index
                    mt_index = i-1
                    # modify mut_df to include mut coor and return new mut_df
                    self.mut_df.at[mt_index, self.mt_row] = self.mt_coor
                    break
                # if i IS NOT between two wt coordinates and it IS NOT the last index in list
                # go to next index in loop
                elif i != index_list[-1]:
                    continue
                # if i IS NOT between two wt coordinates AND it IS the last index in list
                # return an empty df that will cause an exception in run.py                
                else:
                    # here i can add option to look for next start codon (start lost)
                    print(f'{self.tscript_id} mutant junction coordinate is not within coding transcript...skipping')
                    self.mut_df = pd.DataFrame()
                
        # forward direction
        if self.reverse == False:
            # select rows that are possible places to insert mutant coordinate
            revised_df = self.mut_df[self.mut_df[self.wt_row] >= self.wt_coor] 
            # loop over indexes in revised_df
            # ex: [13,14,15] (going to end of df)
            index_list = revised_df.index.tolist()
            for i in index_list:
                # if i IS in between two wt coordinates and i IS NOT the last index in list (i+1 index exists in df)
                if i+1 in index_list and self.mt_coor > revised_df.loc[i, self.wt_row] and self.mt_coor < revised_df.loc[i+1, self.wt_row]:
                    # create mt_index
                    mt_index = i+1
                    # modify mut_df to include mut coor and return new mut_df
                    self.mut_df.at[mt_index, self.mt_row] = self.mt_coor
                    break
                # if i IS NOT between two wt coordinates and it IS NOT the last index in list
                # go to next index in loop
                elif i != index_list[-1]:
                    continue
                # if i IS NOT between two wt coordinates AND it IS the last index in list
                # return an empty df that will cause an exception in run.py
                else:
                    # here i can add option to look for next stop codon (stop lost)
                    print(f'{self.tscript_id} mutant junction coordinate is not within coding transcript...skipping')
                    self.mut_df = pd.DataFrame()    
                
        # exon skip
        # starting knowing both coors are in coding region
        if self.anchor == 'NDA':
            start_index = self.mut_df[self.mut_df[self.wt_row] == self.wt_coor].index.item()
            stop_index = self.mut_df[self.mut_df[self.mt_row] == self.mt_coor].index.item()
            # delete any rows in between
            # [1:] to remove 1st index from index_list (because of range function)
            index_list = list(range(stop_index, start_index))[1:]
            self.mut_df = self.mut_df.drop(index_list)
            
        # testing
        #self.mut_df.to_csv(f'{self.test_dir}/mut_dataframe.tsv', sep='\t', index=False)
        
        return self.mut_df


    def findProtein(self, dataframe):
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

    def createFasta(self, wt_protein, mt_protein):
        #testing
        #file = f'{self.test_dir}/protein.fasta'
        #os.mkdir(f'{self.working_dir}/results')
        file = f'{self.output_dir}/assembled_transcripts.fa'
        write_str = f'>WT.{self.gene_name}.{self.junction_name}.{self.tscript_id}\n{wt_protein}\n>MUT.{self.gene_name}.{self.junction_name}.{self.tscript_id}\n{mt_protein}\n'
        if os.path.exists(file):
            dup_content = re.search(write_str, open(file, "r").read())
            if dup_content == None:
                with open(file, "a") as f:
                    f.write(write_str)
        else:
            with open(file, "a") as f:
                f.write(write_str)


def define_parser():
        parser = argparse.ArgumentParser(
            description="Use Regtools junctions to predict mutant transcript",
        )
        parser.add_argument('fasta', help='path to fasta file', type=str)
        parser.add_argument('transcript', help='Transcript ID', type=str)
        parser.add_argument('chrom', help='Chromosome', type=int)
        parser.add_argument('coordinates', nargs=2, help='Tumor-specific alternative splice junction: start end')
        parser.add_argument('anchor', help='Regtools anchor type', type=str, choices=['A', 'D', 'NDA'])
        parser.add_argument('strand', help='Strand of transcript', type=int)
        parser.add_argument('gene', help='Gene Name', type=str)
        parser.add_argument('output_dir', help='output path', type=str)
        
        return parser

def main(args_input = sys.argv[1:]):
        parser = define_parser()
        args = parser.parse_args(args_input)

        j = JunctionToFasta(args.fasta, args.transcript, args.chrom, args.coordinates, args.anchor, args.strand, args.gene)
        wt = j.createWtDataframe()
        mut = j.createMtDataframe()
        wt_aa = j.findProtein(wt)
        mut_aa = j.findProtein(mut)
        j.createFasta(wt_aa, mut_aa)


if __name__ == '__main__':
    main()
    #j = JunctionToFasta(fasta_path='/Users/mrichters/Documents/ref_fasta/GRCh38.d1.vd1.fa', #tscript_id='ENST00000401026', chrom=1, junction_name='JUNC00000015', junction_coors=[93221940, 93232426], anchor='A', strand=1, gene_name='CCDC18', output_dir='/Users/mrichters/Desktop/Alt_Splicing/HCC1395/results')
    #wt = j.createWtDataframe()
    #mut = j.createMtDataframe()


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
