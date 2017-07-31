import sys, os
import math

#python input_parsing_vaccine_design.py 25 final.tsv
#python input_parsing_vaccine_design.py 25 /gscmnt/gc3018/cancer-genomics/medseq/tmp/mneveau/pvacSeq/epitope_output_12/MHC_Class_I/H_XF-PICI012.final.tsv /gscmnt/gc3018/cancer-genomics/medseq/tmp/mneveau/pvacSeq/H_XF-PICI012.annotated.vcf .

(script, n_mer, input_tsv, input_vcf, output_dir) = sys.argv

output_f = os.path.join(output_dir, "vaccine_design_input.fa")

#get info from pvacseq output final tsv
with open(input_tsv, 'r') as input_f:
    next(input_f)
    mut_IDs, positions, mt_epitope_seqs, wt_epitope_seqs, transcript_IDs = [], [], [], [], []
    for line in input_f:
        fields = line.split("\t")
        mut_type, mutation, pos, gene_name = fields[7], fields[8], fields[9], fields[10]
        mt_epitope_seq, wt_epitope_seq = fields[15], fields[16]
        mutation = mutation.split("/")
        #if position presented as a range, use higher end of range
        if "-" in pos:
            pos = pos.split("-")
            pos = pos[1]
        old_AA, new_AA = mutation[0], mutation[1]
        if mut_type == "missense": 
            mut_ID = "MT." + gene_name + "."  + old_AA + pos + new_AA
        elif mut_type == "FS":
            ##old_AA, new_AA = mutation[0], mt_epitope_seq[0]
            mut_ID = "MT." + gene_name + "." + old_AA + pos + "fs"
        mut_IDs.append(mut_ID)
        ##will probably need for FS
        positions.append(pos)
        mt_epitope_seqs.append(mt_epitope_seq)
        wt_epitope_seqs.append(wt_epitope_seq)
        transcript_IDs.append(fields[5])
input_f.close()

#get info from initial pvacseq input vcf
with open(input_vcf, 'r') as input_f:
    transcripts_dict = {}
    for line in input_f:
        attributes = []
        if line[0] != "#":
            fields = line.split("\t")
            info = fields[7]
            info = info.split("|")
            transcript_ID, downstr_seq, len_change, full_seq = info[6], info[23], info[24], info[25]
            #attributes[0] = full seq, attributes[1] = downstream protein, attributes[2] = protein length change
            attributes.append(full_seq)
            attributes.append(downstr_seq)
            attributes.append(len_change)
            transcripts_dict[transcript_ID] = attributes
    print(transcripts_dict)
input_f.close()

with open(output_f, 'w') as out_f:
    sub_seq = ""
    full_seq = ""
    n_mer = int(n_mer)
    for i in range(len(transcript_IDs)):
        full_seq = (transcripts_dict[transcript_IDs[i]])[0]
        print(transcript_IDs[i])
        print(full_seq)
        #get enough peptides to make n_mer
        beginning = full_seq.find(wt_epitope_seqs[i])
        print(wt_epitope_seqs[i])
        if beginning == -1:
            #deal with frameshift idosyncracies here:
            #must take peptides from front of epitope bc nothing left after
            print("No WT epitope sequence listed")
            downstr_seq, len_change = transcripts_dict[transcript_IDs[i]][1], transcripts_dict[transcript_IDs[i]][2]
            len_needed = n_mer - len(mt_epitope_seqs[i]) 
            start = int(positions[i]) - len_needed
            end = int(positions[i])
            sub_seq = full_seq[start:end] + mt_epitope_seqs[i]
            
            #we want the mt_eptiope sequence + whatever comes before that in the downstr_seq + the remainder needed from the very end of the full_seq

            #whether or not we include the amino acid at the very end (ex R/X vs L/LX) may depend on that 
            
        #find wt sequence and replace with mt seq
        else:
            sub_seq = full_seq.replace(wt_epitope_seqs[i], mt_epitope_seqs[i])
            epitope_len = len(mt_epitope_seqs[i])
            end = beginning + epitope_len
            #if eptitope sequence is too close to the beginning or end to get the
            #right amount of flanking peptides, get appropriate length from solely
            #ahead or behind
            len_needed = n_mer - epitope_len
            print("length needed: " + str(len_needed))
            if len_needed % 2 != 0:
                front = int(beginning - len_needed / 2)
                back = int(end + len_needed / 2) + 1
            else:
                front = beginning - len_needed / 2
                back = end + len_needed / 2
            print(beginning, end, front, back) 
            if front < 0:
                print(1)
                sub_seq = full_seq[beginning:(beginning + n_mer)]
            elif back > len(full_seq):
                print(2)
                sub_seq = full_seq[(end - n_mer):end]
            else:
                print(3)
                print(front)
                print(back)
                sub_seq = full_seq[front:back]
        print("subseq: " + sub_seq)
        print("length: " + str(len(sub_seq)))
        out_f.write(">" + mut_IDs[i] + "\n")
        out_f.write(sub_seq + "\n")
out_f.close()



