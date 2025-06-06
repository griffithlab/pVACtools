import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import tempfile
import re

# Function to break the pepetides ID on the . to extract gene and AA information
def extract_info(value):
    parts = value.split('.')
    result = '.'.join([parts[2], parts[3], parts[4]])
    return result


# Function to rearrange string so that G518D looks like 518G/D
def rearrange_string(s):
    match = re.match(r'([A-Za-z]+)([\d-]+)([A-Za-z]*)', s)
    if match:
        letters_before = match.group(1)
        numbers = match.group(2)
        letters_after = match.group(3)
                
        return f"{numbers}{letters_before}/{letters_after}"
    else:
        return s


# Function to calculate molecular weight---------------------------------------
def calculate_molecular_weight(peptide):
    analyzed_seq = ProteinAnalysis(peptide)
    return analyzed_seq.molecular_weight()


# Function to make id column unique -------------------------------------------
def make_column_unique(df, column_name):
    seen_values = set()
    new_values = []

    for value in df[column_name]:
        if value in seen_values:
            suffix = 1
            while f"{value}.{suffix}" in seen_values:
                suffix += 1
            unique_value = f"{value}.{suffix}"
        else:
            unique_value = value

        seen_values.add(unique_value)
        new_values.append(unique_value)

    df[column_name] = new_values
    return df


def main(peptides_path, classI_path, classII_path, sample_name, output_file, output_path, all_epitopes_flag):
    # Creating the Peptides 51mer Sheet -----------------------------------
    peptides = pd.read_csv(peptides_path, sep="\t")
    peptides =  peptides.drop(['cterm_7mer_gravy_score', 'cysteine_count', 'n_terminal_asparagine', 'asparagine_proline_bond_count', 
                                 'difficult_n_terminal_residue', 'c_terminal_cysteine', 'c_terminal_proline', 'max_7mer_gravy_score'], axis=1)
    peptides["RESTRICTING HLA ALLELE"] = " "

    peptides["CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE MW (CLIENT)"] = peptides["peptide_sequence"].apply(calculate_molecular_weight)

    peptides = peptides.rename(columns={"id":"ID", "peptide_sequence":"CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES"})
    peptides["Comments"] = " "
    peptides["CANDIDATE NEOANTIGEN"] = peptides["ID"].apply(lambda x: '.'.join(x.split('.')[:3]))
    peptides["CANDIDATE NEOANTIGEN"] = sample_name + "." + peptides["CANDIDATE NEOANTIGEN"]

    peptides = peptides[["ID", "CANDIDATE NEOANTIGEN", "CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES", 
                           "RESTRICTING HLA ALLELE", "CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE MW (CLIENT)", "Comments"]]
    
    # Add the Restricting HLA Alles from Class I and Class II
    # create a dataframe that contains the classI and classII pepetide sequence
    # Create a universal ID by editing the peptide 51mer ID
    peptides.rename(columns={'ID': 'full ID'}, inplace=True)
    peptides['51mer ID'] = peptides['full ID']
    peptides['51mer ID'] = peptides['51mer ID'].apply(lambda x: '.'.join(x.split('.')[1:]))  # Removes the 'MT' from the beginning of ID column
    peptides['51mer ID'] = peptides['51mer ID'].apply(lambda x: '.'.join(x.split('.')[1:]))  # Remives the MT index from the ID column
    
    def modify_id(x):
        parts = x.split('.')
        
        def is_valid_char(char):
            return char.isdigit() or char == '-' or (char in ['.', '/'] and parts[-1][parts[-1].index(char)+1].isdigit())
        
        if not 'missense' in parts:
            if 'FS' in parts:
                # MT.239.KIF7.ENST00000394412.8.FS.405-410CGCGCACTCGGCGCCCAG/C
                # KIF7.ENST00000394412.FS405-410
                # If 'FS' is present, remove non-digit characters after the last period
                last_part = ''.join(filter(is_valid_char, parts[-1]))
                modified_id = '.'.join(parts[:-1]) + last_part
            elif 'inframe_del' in parts:
                # MT.328.DNAAF3.ENST00000391720.8.inframe_del.585-589KTGV*/R
                # DNAAF3.ENST00000391720.8.KTGV*538-542R
                last_part = parts[-1]
                position = re.match(r'[^a-zA-Z]*', last_part).group()
                aa = re.search(r'[^0-9-]+.*', last_part).group()
                aa_parts = aa.split('/')
                modified_id = '.'.join(parts[:3])
                aa_change = aa_parts[0] + position + aa_parts[1]
                modified_id = modified_id + "." + aa_change
            elif 'inframe_ins' in parts:
                # MT.249.ZCCHC14.ENST00000268616.9.inframe_ins.769H/HH
                # ZCCHC14.ENST00000268616.9.H769HH
                modified_id = '.'.join(parts[:3] + parts[4:])
            else:
                print("Non missense candidate not accounted for!")
                print(parts)     
        else:
            # If 'missense' or other labels are present, remove them
            # MT.244.ZP2.ENST00000574002.1.missense.84D/Y
            # ZP2.ENST00000574002.1.84D/Y
            modified_id = '.'.join(parts[:3] + parts[4:])
        
        return modified_id
    
    

    #peptides['51mer ID'] = peptides['51mer ID'].apply(lambda x: '.'.join(x.split('.')[:3]) + '.' + '.'.join(x.split('.')[4:])) # removes the variant label
    peptides['51mer ID'] = peptides['51mer ID'].apply(modify_id)

    classI = pd.read_csv(classI_path, sep="\t")
    classII = pd.read_csv(classII_path, sep="\t")
    if all_epitopes_flag:
        classI.rename(columns = {"MT Epitope Seq":"Best Peptide Class I", "HLA Allele":"Class I Allele", 
                                 "Median MT IC50 Score":"Class I IC50 MT", "Median MT Percentile":"Class I %ile MT", 
                                 "Transcript":"Class I Best Transcript"}, inplace=True)
        classII.rename(columns = {"MT Epitope Seq":"Best Peptide Class II", "HLA Allele":"Class II Allele", 
                                  "Median MT IC50 Score":"Class II IC50 MT", "Median MT Percentile":"Class II %ile MT", 
                                  "Transcript":"Class II Best Transcript"}, inplace=True)
        
        classI['position AA Change'] = classI['Index'].split('.')[5]
        classI['51mer ID'] = classI['Gene Name'] + '.' + classI['Class I Best Transcript'] + '.' + classI['position AA Change'] 
        class_sequences = pd.merge(classI[['Index', 'Best Peptide Class I', '51mer ID', 'Pos', 'AA Change', 'Class I Allele', "Class I IC50 MT", "Class I %ile MT", "Class I Best Transcript"]], 
                                   classII[['Index', 'Best Peptide Class II', 'Class II Allele', "Class II IC50 MT", "Class II %ile MT", "Class II Best Transcript"]], on='ID', how='left')
        class_sequences = class_sequences.drop(columns=['Index'])

    else:
        classI.rename(columns = {"Best Peptide":"Best Peptide Class I", "Allele":"Class I Allele", 
                                 "IC50 MT":"Class I IC50 MT", "%ile MT":"Class I %ile MT", 
                                 "Best Transcript":"Class I Best Transcript"}, inplace=True)
        classII.rename(columns = {"Best Peptide":"Best Peptide Class II", "Allele":"Class II Allele", 
                                  "IC50 MT":"Class II IC50 MT", "%ile MT":"Class II %ile MT", 
                                  "Best Transcript":"Class II Best Transcript"}, inplace=True)
    
        def rearrange_string(s):
            match = re.match(r'([A-Za-z]+)([\d-]+)([A-Za-z]+)', s)
            if match:
                letters_before = match.group(1)
                numbers = match.group(2)
                letters_after = match.group(3)
                    
                return f"{numbers}{letters_before}/{letters_after}"
                    # Just use the postion for the key to avoid FS problem
                #return f"{numbers}"
            else:
                return s
        
    
        classI['position AA Change'] = classI['AA Change'].apply(rearrange_string)
        classI['51mer ID'] = classI['Gene'] + '.' + classI['Class I Best Transcript'] + '.' + classI['position AA Change'] 
        class_sequences = pd.merge(classI[['ID', 'Best Peptide Class I', '51mer ID', 'Pos', 'AA Change', 'Class I Allele', "Class I IC50 MT", "Class I %ile MT", "Class I Best Transcript"]], 
                                   classII[['ID', 'Best Peptide Class II', 'Class II Allele', "Class II IC50 MT", "Class II %ile MT", "Class II Best Transcript"]], on='ID', how='left')
        class_sequences = class_sequences.drop(columns=['ID'])

    
    merged_peptide_51mer = pd.merge(peptides, class_sequences, on='51mer ID', how='left')
    
    merged_peptide_51mer['sorting id'] = merged_peptide_51mer['full ID'].apply(extract_info) # creating a ID to sort reviewed canidates by the order of the 51mer
    merged_peptide_51mer = make_column_unique(merged_peptide_51mer, 'sorting id') # make sure every sorting id is unique

    # Dropping the sorting column -------------------------
    merged_peptide_51mer = merged_peptide_51mer.drop(columns=['sorting id'])

    with tempfile.NamedTemporaryFile(suffix="_Peptides_51-mer.xlsx", delete=False) as tmp:
        Peptide_file_name = tmp.name
    merged_peptide_51mer.to_excel(Peptide_file_name, index=False)
    
    return Peptide_file_name


if __name__ == "__main__":
    main()
