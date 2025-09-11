import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import tempfile
import vcfpy
import os
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


# Function that reads a vcf and returns the variants that contain "PASS" in "FILTER"
def load_pass_variants(vcf_path):
    variants = set()
    reader = vcfpy.Reader.from_path(vcf_path)
    for record in reader:
        if record.FILTER and record.FILTER[0] != "PASS":
            continue
        chrom = record.CHROM
        ref = record.REF
        for alt in record.ALT:
            alt_str = alt.value
            start = record.affected_start
            end = record.affected_end
            variants.add((chrom, start, end, ref, alt_str))
    return variants


# Function that fills the "Variant Called in External VCF" column based on the presence
# of the classI_tsv variant in the provided VCF along with a value of "PASS" in "FILTER"
def fill_variant_called_column(df, input_vcf, external_vcf):
    input_variants = load_pass_variants(input_vcf)
    external_variants = load_pass_variants(external_vcf)

    called = []
    for _, row in df.iterrows():
        parts = row["ID"].split("-")
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        ref = parts[3]
        alt = parts[4]
        variant_key = (chrom, start, end, ref, alt)

        if variant_key in input_variants and variant_key in external_variants:
            called.append(True)
        else:
            called.append(False)
    df["Variant Called in External VCF"] = called


def main(peptides_path, classI_path, classII_path, input_vcf, external_vcf, sample_name,
        allowed_evaluations, output_file_prefix, output_path):
    reviewed_candidates = pd.read_csv(classI_path, sep="\t")
    reviewed_candidates = reviewed_candidates[reviewed_candidates["Evaluation"].isin(allowed_evaluations)]

    reviewed_candidates = reviewed_candidates.rename(columns={'Comments':'pVAC Review Comments'})
    reviewed_candidates["IGV Review Comments"] = ""

    if external_vcf:
        fill_variant_called_column(reviewed_candidates, input_vcf, external_vcf)

    # create sorting ID that is gene and transcript to sort in the same order as peptide
    reviewed_candidates['sorting id'] = reviewed_candidates['Gene']  + '.' + reviewed_candidates['Best Transcript']
    # make sure the sorting id column is unique
    reviewed_candidates = make_column_unique(reviewed_candidates, 'sorting id')

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
    peptides['51mer ID'] = peptides['51mer ID'].apply(lambda x: x.split('.', 1)[1])  # Removes the 'MT' from the beginning of ID column
    peptides['51mer ID'] = peptides['51mer ID'].apply(lambda x: x.split('.', 1)[1])  # Remives the MT index from the ID column

    def modify_id(original_id):
        match = re.match(r'([\w\-]+)\.(ENS[0-9A-Z]+\d+(?:\.\d+)?)\.(\w+)\.(.+)$', original_id)
        if not match:
            raise Exception(f"Unrecognized ID format: {original_id}")
            return original_id

        gene, transcript, variant_type, last_part = match.groups()

        if variant_type == "missense" or variant_type == "inframe_ins":
            # Removes the variant type from the ID
            # MT.244.ZP2.ENST00000574002.1.missense.84D/Y
            # ZP2.ENST00000574002.1.84D/Y
            modified_id = f"{gene}.{transcript}.{last_part}"

        elif variant_type == "FS":
            # MT.239.KIF7.ENST00000394412.8.FS.405-410CGCGCACTCGGCGCCCAG/C
            # KIF7.ENST00000394412.8.FS405-410
            fs_pos = re.sub(r'[^\d\-]', '', last_part)
            modified_id = f"{gene}.{transcript}.FS{fs_pos}"

        elif variant_type == "inframe_del":
            # MT.328.DNAAF3.ENST00000391720.8.inframe_del.585-589KTGV*/R
            # DNAAF3.ENST00000391720.8.KTGV*538-542R
            m = re.match(r'(\d+-\d+)([A-Z\*/]+)\/([A-Z\*/]+)', last_part)
            if m:
                pos, ref, alt = m.groups()
                modified_id = f"{gene}.{transcript}.{ref}{pos}{alt}"
            else:
                print(f"Unrecognized inframe_del format: {original_id}")
                modified_id = original_id

        else:
            print(f"Non-missense variant not handled: {original_id}")
            modified_id = original_id

        return modified_id


    #peptides['51mer ID'] = peptides['51mer ID'].apply(lambda x: '.'.join(x.split('.')[:3]) + '.' + '.'.join(x.split('.')[4:])) # removes the variant label
    peptides['51mer ID'] = peptides['51mer ID'].apply(modify_id)

    classI = pd.read_csv(classI_path, sep="\t")
    classII = pd.read_csv(classII_path, sep="\t")

    classI.rename(columns = {"Best Peptide":"Best Peptide Class I", "Allele":"Class I Allele",
                                "IC50 MT":"Class I IC50 MT", "%ile MT":"Class I %ile MT",
                                "Best Transcript":"Class I Best Transcript"}, inplace=True)
    classII.rename(columns = {"Best Peptide":"Best Peptide Class II", "Allele":"Class II Allele",
                                "IC50 MT":"Class II IC50 MT", "%ile MT":"Class II %ile MT",
                                "Best Transcript":"Class II Best Transcript"}, inplace=True)

    def rearrange_string(s):
        # Handle insertions like -287-288L → 287-288-/L
        insertion_match = re.match(r'-?(\d+-\d+)([A-Za-z]+)', s)
        if insertion_match:
            position = insertion_match.group(1)
            inserted = insertion_match.group(2)
            return f"{position}-/{inserted}"

        # Handle substitutions like E510Q → 510E/Q
        substitution_match = re.match(r'([A-Za-z]+)(\d+)([A-Za-z]+)', s)
        if substitution_match:
            letters_before = substitution_match.group(1)
            numbers = substitution_match.group(2)
            letters_after = substitution_match.group(3)
            return f"{numbers}{letters_before}/{letters_after}"

        # Return unchanged if no pattern matched
        return s


    classI['position AA Change'] = classI['AA Change'].apply(rearrange_string)
    classI['51mer ID'] = classI['Gene'] + '.' + classI['Class I Best Transcript'] + '.' + classI['position AA Change'] 
    class_sequences = pd.merge(classI[['ID', 'Best Peptide Class I', '51mer ID', 'Pos', 'AA Change', 'Class I Allele', "Class I IC50 MT", "Class I %ile MT", "Class I Best Transcript"]], 
                                classII[['ID', 'Best Peptide Class II', 'Class II Allele', "Class II IC50 MT", "Class II %ile MT", "Class II Best Transcript"]], on='ID', how='left')
    class_sequences = class_sequences.drop(columns=['ID'])

    merged_peptide_51mer = pd.merge(peptides, class_sequences, on='51mer ID', how='left')

    merged_peptide_51mer['sorting id'] = merged_peptide_51mer['full ID'].apply(extract_info) # creating a ID to sort reviewed canidates by the order of the 51mer
    merged_peptide_51mer = make_column_unique(merged_peptide_51mer, 'sorting id') # make sure every sorting id is unique

    # Sorting Candidates sheet to be in the same order as Peptides Sheet -------------------------
    reviewed_candidates = reviewed_candidates.set_index('sorting id')
    reviewed_candidates = reviewed_candidates.reindex(index=merged_peptide_51mer['sorting id'])
    reviewed_candidates = reviewed_candidates.reset_index()

    # Dropping the sorting column -------------------------
    reviewed_candidates = reviewed_candidates.drop(columns=['sorting id'])
    merged_peptide_51mer = merged_peptide_51mer.drop(columns=['sorting id'])

    with tempfile.NamedTemporaryFile(suffix="_Peptides_51-mer.xlsx", delete=False) as tmp:
        Peptide_file_name = tmp.name
    merged_peptide_51mer.to_excel(Peptide_file_name, index=False)

    reviewed_candidates_file_name = f"{output_file_prefix}_{sample_name}.Annotated.Neoantigen_Candidates.xlsx"
    reviewed_candidates_file_path = os.path.join(output_path, reviewed_candidates_file_name)
    reviewed_candidates.to_excel(reviewed_candidates_file_path, index=False)

    return Peptide_file_name


if __name__ == "__main__":
    main()
