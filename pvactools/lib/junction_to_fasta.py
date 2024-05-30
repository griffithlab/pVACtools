import os
import re
import pandas as pd
import pyfaidx
from Bio.Seq import Seq

class JunctionToFasta():
    def __init__(self, **kwargs):
        self.fasta_path     = kwargs['fasta_path']
        self.tscript_id     = kwargs['tscript_id']
        self.chrom          = kwargs['chrom']
        self.gtf_df         = kwargs['gtf_df']
        self.junction_name  = kwargs['junction_name']
        self.junction_coors = kwargs['junction_coors']
        self.junction_df    = kwargs['junction_df']
        self.fasta_index    = kwargs['fasta_index']
        self.variant_info   = kwargs['variant_info']
        self.anchor         = kwargs['anchor']
        self.strand         = kwargs['strand']
        self.gene_name      = kwargs['gene_name']
        self.output_file    = kwargs['output_file']
        self.output_dir     = kwargs['output_dir']
        self.sample_name    = kwargs['sample_name']
        self.vcf_file       = kwargs['vcf']
        if (self.anchor == 'A' and self.strand == 1) or (self.anchor == 'D' and self.strand == -1) or (self.anchor == 'NDA'):
            self.wt_coor  = int(self.junction_coors[1])
            self.alt_coor = int(self.junction_coors[0])
            self.wt_row   = "cds_start"
            self.alt_row  = "cds_stop"
            self.reverse  = True
        elif (self.anchor == 'D' and self.strand == 1) or (self.anchor == 'A' and self.strand == -1):
            self.wt_coor  = int(self.junction_coors[0])
            self.alt_coor = int(self.junction_coors[1])
            self.wt_row   = "cds_stop"
            self.alt_row  = "cds_start"
            self.reverse  = False
        self.personal_fasta = pyfaidx.FastaVariant(self.fasta_path, self.vcf_file, sample=self.sample_name)

    def load_gtf_data(self):
        # subset df by transcript_id, get coding coordinates
        # sort by cds_start - if on reverse strand, will be reverse comp. in aa sequence
        coding_df = self.gtf_df[(self.gtf_df['transcript_id'] == self.tscript_id) & (self.gtf_df['feature'] == 'CDS')][['cds_chrom', 'cds_start', 'cds_stop']].sort_values('cds_start').reset_index(drop=True)

        return coding_df

    def create_wt_df(self):
        # load in wt transcript df
        self.wt_df = self.load_gtf_data()
        # if anchor is D or A, make sure the wt coordinate is inside in the coding region of transcript
        # (alt coordinate won't be present because ensembl only lists the wt coordinates) 
        if self.anchor in ['D', 'A']:
            # look for ref coor index
            index = [index for index,value in self.wt_df[self.wt_row].items() if value == self.wt_coor]
            # value not found in df (Exception)
            if not index:
                print(f'{self.tscript_id} WT coordinate is not within coding transcript...Skipping')
                print(f'Missing coordinate ({self.anchor}): {self.wt_coor}')
                self.wt_df = pd.DataFrame()
            # else:
            #     print(f'{self.anchor} WT: {self.wt_coor} {index[0]}')
        # if exon skip, check that both coordinates are inside the coding region of transcript
        elif self.anchor == 'NDA':                    
            # check for presence of both coordinates 
            index_wt = [index for index,value in self.wt_df[self.wt_row].items() if value == self.wt_coor]
            index_alt = [index for index,value in self.wt_df[self.alt_row].items() if value == self.alt_coor]
            if index_wt and index_alt:
                # print(f'NDA both coors present: {self.wt_coor}, {index_wt[0]}, {self.alt_coor}, {index_alt[0]}')
                pass
            else:    
                if not index_wt and not index_alt:
                    print(f'{self.fasta_index} Exon skip: both junction coordinates outside of coding transcript...Skipping')
                    print(f'Missing coordinates: {self.wt_coor} {self.alt_coor}')
                elif (not index_wt and index_alt) or (index_wt and not index_alt):
                    print(f'{self.fasta_index} Exon skip: one junction coordinate outside of coding transcript...Skipping')
                    if index_wt:
                        print(f'Missing coordinate: {self.alt_coor}')
                    elif index_alt:
                        print(f'Missing coordinate: {self.wt_coor}')
                self.wt_df = pd.DataFrame()
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
                    #print(f'{self.anchor} ALT: {self.alt_coor} {alt_index}')
                    break
                # if i IS NOT between two wt coordinates and it IS NOT the last index in list
                # go to next index in loop
                elif i != index_list[-1]:
                    continue
                # if i IS NOT between two wt coordinates AND it IS the last index in list
                # return an empty df that will cause an exception in run.py                
                else:
                    # here i can add option to look for next start codon (start lost)
                    print(f'{self.fasta_index} Alternate junction coordinate outside of coding transcript...Skipping')
                    print(f'Missing coordinate: {self.alt_coor}')
                    self.alt_df = pd.DataFrame()
                    continue
                
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
                    #print(f'{self.anchor} ALT: {self.alt_coor} {alt_index}')
                    break
                # if i IS NOT between two wt coordinates and it IS NOT the last index in list
                # go to next index in loop
                elif i != index_list[-1]:
                    continue
                # if i IS NOT between two wt coordinates AND it IS the last index in list
                # return an empty df that will cause an exception in run.py
                else:
                    # here i can add option to look for next stop codon (stop lost)
                    print(f'{self.fasta_index} Alternate junction coordinate outside of coding transcript...Skipping')
                    print(f'Missing coordinate: {self.alt_coor}')
                    self.alt_df = pd.DataFrame()
                    continue    
                
        # exon skip
        # starting knowing both coors are in coding region
        if self.anchor == 'NDA':
            start_index = self.alt_df[self.alt_df[self.wt_row] == self.wt_coor].index.item()
            stop_index = self.alt_df[self.alt_df[self.alt_row] == self.alt_coor].index.item()
            # delete any rows in between
            # [1:] to remove 1st index from index_list (because of range function)
            index_list = list(range(stop_index, start_index))[1:]
            self.alt_df = self.alt_df.drop(index_list)
        return self.alt_df

    def get_aa_sequence(self, dataframe):
        # create coding_coors column for fasta indexing
        dataframe["coding_coors"] = dataframe["cds_start"].astype(str) + "," + dataframe["cds_stop"].astype(str)
        coordinates = dataframe["coding_coors"].tolist()
        # generate AA sequence from coding exon coordinates (pyfaidx)
        # taking a lot longer when using a gzipped personal fasta - why?
        final_seq = ''
        for x in coordinates:
            start = int(x.split(',')[0]); end = int(x.split(',')[1])
            seq = self.personal_fasta.get_seq(self.chrom, start, end)
            final_seq += str(seq)
        # using Seq from Bio.Seq to translate str_seq
        # positive strand
        dna_seq = Seq(final_seq) # creating Seq object
        # negative strand
        if self.strand == -1:
            dna_seq = dna_seq.reverse_complement() # still a Seq object
        # make biopython happy by making dna_seq a multiple of 3 
        # adding Ns to end of sequence if remainder != 0
        remainder = len(dna_seq) % 3
        if remainder == 0:
            frameshift = 'no'
        else:
            frameshift = 'yes'
            if remainder == 1:
                dna_seq += "NN"
            elif remainder == 2:
                dna_seq += "N"
        aa_seq = str(dna_seq.translate(to_stop=True)) # translating from Seq object
        if aa_seq[0] != 'M':
            print(f'{self.tscript_id} does not begin with start codon...Skipping')
            aa_seq = ''

        return aa_seq, frameshift

    def create_sequence_fasta(self, wt_seq, alt_seq):
        write_str = f'>WT.{self.fasta_index}\n{wt_seq}\n>ALT.{self.fasta_index}\n{alt_seq}\n'
        if os.path.exists(self.output_file):
            with open(self.output_file, 'r+') as f:
                dup_content = re.search(write_str, f.read())
                if dup_content == None:
                    f.write(write_str)
        else:
            with open(self.output_file, 'w') as e:
                e.write(write_str)

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
