import pandas as pd
from bs4 import BeautifulSoup

class AminoAcid:

    def __init__(self, nucleotide, bold, color, underline, large, position, open_tag, close_tag):
        self.nucleotide = nucleotide
        self.bold = bold
        self.color = color
        self.underline = underline
        self.large = large
        self.position = position 
        self.open_tag = open_tag
        self.close_tag = close_tag

    def view(self):
        print("Nucleotide: ", self.nucleotide)
        print("Open Tag: ", self.open_tag)
        print("Close Tag: ", self.close_tag)
        print("Bold: ", self.bold)
        print("Color: ",self.color)
        print("Underline: ", self.underline)
        print("Large: ", self.large)


def annotate_every_nucleotide(sequence, classI_peptide, classII_peptide, 
                              classI_ic50, classI_percentile, classII_ic50, classII_percentile, 
                              classI_transcript, classII_transcript, 
                              cIIC50_threshold, cIpercentile_threshold, cIIIC50_threshold, cIIpercent_threshold, probPos):

    peptide_sequence = [] # Create a list to hold all AA of the 51mer

    # Make the sequence a list of AminoAcid objects
    for i in range(len(sequence)):
        new_AA = AminoAcid(sequence[i], False, False, False, False, -1, False, False)
        
        if len(probPos) > 0 and sequence[i] in probPos:
            new_AA.large = True
        
        peptide_sequence.append(new_AA)

    # CLASS I
    positions = []
    # Get the positions in the peptide_sequence where the classI is located
    for i in range(len(peptide_sequence)):
        for j in range(len(classI_peptide)):
            if peptide_sequence[i].nucleotide == classI_peptide[j]:
                positions.append(i)
                i+=1
            else:
                break

        if len(positions) == len(classI_peptide):
            break
        else:
            positions = []
            
    # set those positions to red
    # if median affinity < 1000 nm OR percentile < 2%
    j = 0
    if float(classI_ic50) < cIIC50_threshold or float(classI_percentile) < cIpercentile_threshold:
        for i in range(len(peptide_sequence)):
            if j < len(positions) and i == positions[j]:
                peptide_sequence[i].color = True
                j+=1

    if classI_transcript == classII_transcript:
        # CLASS II         
        positions = []
        for i in range(len(peptide_sequence)):
            for j in range(len(classII_peptide)):
                if peptide_sequence[i].nucleotide == classII_peptide[j]:
                    positions.append(i)
                    i+=1
                else:
                    break

            if len(positions) == len(classII_peptide):
                break
            else:
                positions = []
        # Set class II to bold
        # if percentile < 2% 
        j = 0
        if float(classII_percentile) < cIIpercent_threshold or float(classII_ic50) < cIIIC50_threshold:
            for i in range(len(peptide_sequence)):
                if j < len(positions) and i == positions[j]:
                    peptide_sequence[i].bold = True
                    j+=1
    else:
        print("Note: ClassII transcript different then ClassI. ClassII peptide not bolded.")
    

    return(peptide_sequence)

def set_underline(peptide_sequence, mutant_peptide_pos, row_ID):

    frameshift = False
    classI_position = 0

    # Determine if frameshift mutation by seraching for = '-'
    if ',' in mutant_peptide_pos:
        positions = mutant_peptide_pos.split(",")
        start_position = int(positions[0])
        end_position = int(positions[1])
        frameshift = True
    elif '-' in mutant_peptide_pos: # to account for older versions of pvacseq
        positions = mutant_peptide_pos.split("-")
        start_position = int(positions[0])
        end_position = int(positions[1])
        frameshift = True
    elif mutant_peptide_pos == 'nan': 
        return
    else:
        mutant_peptide_pos = int(float(mutant_peptide_pos))

    if frameshift:
        continue_underline = False
        
        for i in range(len(peptide_sequence)):

            if peptide_sequence[i].color:
                classI_position += 1
            else:
                 classI_position = 0
                 continue_underline = False

            if classI_position == start_position:
                peptide_sequence[i].underline = True
                continue_underline = True
            elif continue_underline:
                peptide_sequence[i].underline = True
            elif classI_position == end_position:
                peptide_sequence[i].underline = True
                continue_underline = False
            i+=1
    else:
        for i in range(len(peptide_sequence)):

            if peptide_sequence[i].color:
                classI_position += 1
            else:
                 classI_position = 0

            if classI_position == int(mutant_peptide_pos):
                peptide_sequence[i].underline = True
            i+=1
        
def set_span_tags(peptide_sequence):

    currently_bold = False
    currently_red = False
    currently_underlined = False
    currently_large = False
    inside_span = False

    for nucleotide in peptide_sequence:

        if currently_bold != nucleotide.bold or currently_red != nucleotide.color or currently_underlined != nucleotide.underline or currently_large != nucleotide.large:
            
            nucleotide.open_tag = True

            if inside_span:
                nucleotide.close_tag = True # only if its inside a span tag
            else:
                nucleotide.close_tag = False


            currently_bold = nucleotide.bold 
            currently_red = nucleotide.color
            currently_underlined = nucleotide.underline
            currently_large = nucleotide.large

            inside_span = True
        
    return(peptide_sequence)

def create_stylized_sequence(peptide_sequence):

    new_string = ''

    for nucleotide in peptide_sequence:

        if nucleotide.open_tag or nucleotide.close_tag:
            if nucleotide.close_tag:
                new_string += '</span>'
                
            if nucleotide.open_tag:
                new_string += '<span style="'
                if nucleotide.bold:
                    new_string += 'font-weight:bold;'
                if nucleotide.color:
                    new_string += 'color:#ff0000;'
                if nucleotide.underline:
                    new_string += 'text-decoration:underline;'
                if nucleotide.large:
                     new_string += 'font-size:105%;'
                new_string += '">'
                new_string += nucleotide.nucleotide
        else:
            new_string += nucleotide.nucleotide
    return(new_string)

def main(peptides_path, sample_name, classI_ic50_score_max, classI_ic50_percentile_max, 
         classII_ic50_score_max, classII_ic50_percentile_max, problematic_position, output_file, output_path):
    # read in classI and class II
    peptides_51mer = pd.read_excel(peptides_path)
 
    # Fill in the Restricting HLA Allele Column
    for index, row in peptides_51mer.iterrows():
        restricting_alleles = ''
        if (float(row['Class I IC50 MT']) < classI_ic50_score_max or float(row['Class I %ile MT']) <  classI_ic50_percentile_max):
            restricting_alleles = row['Class I Allele']
        
        if (float(row['Class I IC50 MT']) < classII_ic50_score_max or float(row['Class I %ile MT']) <  classII_ic50_percentile_max):
            
            if restricting_alleles != '':
                restricting_alleles = restricting_alleles + '/' + row['Class II Allele']
            else:
                 restricting_alleles =  row['Class II Allele']
        
        peptides_51mer.at[index, 'RESTRICTING HLA ALLELE'] = restricting_alleles
                

    # convert peptide 51mer to HTML
    peptides_51mer_html = peptides_51mer.to_html(index=False) # convert to html

    # Creating a BeautifulSoup object and specifying the parser
    peptides_51mer_soup = BeautifulSoup(peptides_51mer_html, 'html.parser')

    for index, row in peptides_51mer.iterrows():
        search_string = row['51mer ID']

        # classII sequence 
        classII_peptide = peptides_51mer.loc[peptides_51mer['51mer ID'] == search_string, 'Best Peptide Class II'].values[0]
        # classI sequence 
        classI_peptide = peptides_51mer.loc[peptides_51mer['51mer ID'] == search_string, 'Best Peptide Class I'].values[0]
        # mutant peptide position
        mutant_peptide_pos = str(peptides_51mer.loc[peptides_51mer['51mer ID'] == search_string, 'Pos'].values[0])
        # classI IC50
        classI_ic50 = peptides_51mer.loc[peptides_51mer['51mer ID'] == search_string, 'Class I IC50 MT'].values[0]
        # classI percentile
        classI_percentile = peptides_51mer.loc[peptides_51mer['51mer ID'] == search_string, 'Class I %ile MT'].values[0]
        # classII IC50
        classII_ic50 = peptides_51mer.loc[peptides_51mer['51mer ID'] == search_string, 'Class II IC50 MT'].values[0]
        # classII percentile
        classII_percentile = peptides_51mer.loc[peptides_51mer['51mer ID'] == search_string, 'Class II %ile MT'].values[0]
        # classI transcript
        classI_transcript = peptides_51mer.loc[peptides_51mer['51mer ID'] == search_string, 'Class I Best Transcript'].values[0]
        # classI transcript
        classII_transcript = peptides_51mer.loc[peptides_51mer['51mer ID'] == search_string, 'Class II Best Transcript'].values[0]


        # Find the tag containing the search string
        tag_with_search_string = peptides_51mer_soup.find('td', string=search_string)

        if tag_with_search_string and isinstance(classII_peptide, str):

            # Find the parent <tr> tag of the tag containing the search string
            parent_tr = tag_with_search_string.find_parent('tr')    
            # Find the next two <td> tags
            next_td_tags = parent_tr.findChildren('td', limit=3)
            
            sequence = next_td_tags[2].get_text()

            # make sequence the list of objects
            peptide_sequence = annotate_every_nucleotide(sequence, classI_peptide, classII_peptide, 
                                                         classI_ic50, classI_percentile, classII_ic50, classII_percentile,
                                                         classI_transcript, classII_transcript,
                                                         classI_ic50_score_max, classI_ic50_percentile_max,
                                                         classII_ic50_score_max, classII_ic50_percentile_max, problematic_position)

            # actually lets break class I and classII into two steps and handle the mutated nucleotide in class I function
            # it should be basically like at that position in the class I set 
            
            set_underline(peptide_sequence, mutant_peptide_pos, row['51mer ID'])

            set_span_tags(peptide_sequence) # pass by reference
            
            print(row['51mer ID'])
            new_string = create_stylized_sequence(peptide_sequence)

            next_td_tags[2].string = new_string

            # Remove the tag_with_search_string from the BeautifulSoup tree
            tag_with_search_string.decompose()


        else:
            print("\nNOT FOUND: ", search_string)
            print("Mutant Peptide Position: ", mutant_peptide_pos)
            print("ClassI: ", classI_peptide)
            print("ClassII: ", classII_peptide, "\n")


        tag_with_search_string = peptides_51mer_soup.select_one('th:-soup-contains("51mer ID")')
        if tag_with_search_string:
            tag_with_search_string.decompose()
        # Now 'soup' contains the modified HTML with the tag removed
        modified_html = peptides_51mer_soup.prettify(formatter=None)

        print()

    html_file_name = f"{output_path}{output_file}_{sample_name}.Colored_Peptides.html" 

    with open(html_file_name, "w", encoding = 'utf-8') as file:
        file.write(modified_html)

if __name__ == "__main__":
    main()
