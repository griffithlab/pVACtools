import pandas as pd
import os
import xlsxwriter


class AminoAcid:
    def __init__(self, nucleotide, bold, color, underline, large, position):
        self.nucleotide = nucleotide
        self.bold = bold
        self.color = color
        self.underline = underline
        self.large = large
        self.position = position

    def view(self):
        print("Nucleotide: ", self.nucleotide)
        print("Bold: ", self.bold)
        print("Color: ", self.color)
        print("Underline: ", self.underline)
        print("Large: ", self.large)


def annotate_every_nucleotide(
    sequence,
    classI_peptide,
    classII_peptide,
    classI_ic50,
    classI_percentile,
    classII_ic50,
    classII_percentile,
    classI_transcript,
    classII_transcript,
    cIIC50_threshold,
    cIpercentile_threshold,
    cIIIC50_threshold,
    cIIpercent_threshold,
    probPos,
):

    peptide_sequence = []

    for i in range(len(sequence)):
        large = sequence[i] in probPos if probPos else False
        new_AA = AminoAcid(sequence[i], False, False, False, large, -1)
        peptide_sequence.append(new_AA)

    sequence_str = "".join(aa.nucleotide for aa in peptide_sequence)
    start_index = sequence_str.find(classI_peptide)

    if start_index != -1:
        positions = list(range(start_index, start_index + len(classI_peptide)))
    else:
        positions = []

    if (
        float(classI_ic50) < cIIC50_threshold
        or float(classI_percentile) < cIpercentile_threshold
    ):
        for position in positions:
            peptide_sequence[position].color = True

    if classI_transcript == classII_transcript:
        start_index = sequence_str.find(classII_peptide)
        if start_index != -1:
            positions = list(range(start_index, start_index + len(classII_peptide)))
        else:
            positions = []
        if (
            float(classII_percentile) < cIIpercent_threshold
            or float(classII_ic50) < cIIIC50_threshold
        ):
            for position in positions:
                peptide_sequence[position].bold = True
    else:
        print(
            "Note: ClassII transcript different than ClassI. ClassII peptide not bolded."
        )

    return peptide_sequence


def set_underline(peptide_sequence, mutant_peptide_pos, full_row_ID):
    frameshift = False
    classI_position = 0

    if ".FS." in full_row_ID:
        frameshift = True
    elif mutant_peptide_pos == "nan":
        return
    else:
        mutant_peptide_pos = int(float(mutant_peptide_pos))

    if frameshift:
        # Find last mutated (colored) position
        last_mut_idx = -1
        for i in reversed(range(len(peptide_sequence))):
            if peptide_sequence[i].color:
                last_mut_idx = i
                break

        if last_mut_idx == -1:
            # No colored residues, nothing to underline
            pass
        else:
            # Walk backward to find contiguous colored positions
            start_idx = last_mut_idx
            for i in reversed(range(0, last_mut_idx)):
                if peptide_sequence[i].color:
                    start_idx = i
                else:
                    break  # stop at first non-colored amino acid

            # Walk forward to find contiguous colored positions
            end_idx = last_mut_idx
            for i in range(last_mut_idx + 1, len(peptide_sequence)):
                if peptide_sequence[i].color:
                    end_idx = i
                else:
                    break  # stop at first non-colored amino acid

            # Underline the identified range
            for i in range(start_idx, end_idx + 1):
                peptide_sequence[i].underline = True
    else:
        for i in range(len(peptide_sequence)):
            if peptide_sequence[i].color:
                classI_position += 1
            else:
                classI_position = 0
            if classI_position == int(mutant_peptide_pos):
                peptide_sequence[i].underline = True


def generate_formatted_excel(peptides_df, output_path, output_file, sample_name):
    file_name = f"{output_file}_{sample_name}.Colored_Peptides.xlsx"
    file_path = os.path.join(output_path, file_name)
    workbook = xlsxwriter.Workbook(file_path)
    worksheet = workbook.add_worksheet()

    bold_fmt = workbook.add_format({"bold": True})
    red_fmt = workbook.add_format({"font_color": "red"})
    bold_red_fmt = workbook.add_format({"bold": True, "font_color": "red"})
    red_underline_fmt = workbook.add_format({"font_color": "red", "underline": True})
    bold_red_underline_fmt = workbook.add_format(
        {"bold": True, "font_color": "red", "underline": True}
    )

    visible_columns = [col for col in peptides_df.columns if col != "Stylized Sequence"]
    for col, header in enumerate(visible_columns):
        worksheet.write(0, col, header)

    for row_idx, row in peptides_df.iterrows():
        for col_idx, header in enumerate(visible_columns):
            value = row[header]
            if (
                header
                == "CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES"
            ):
                peptide_sequence = row["Stylized Sequence"]
                segments = []
                for aa in peptide_sequence:
                    fmt = None
                    if aa.bold and aa.color and aa.underline:
                        fmt = bold_red_underline_fmt
                    elif aa.color and aa.underline:
                        fmt = red_underline_fmt
                    elif aa.bold and aa.color:
                        fmt = bold_red_fmt
                    elif aa.bold:
                        fmt = bold_fmt
                    elif aa.color:
                        fmt = red_fmt

                    if fmt:
                        segments.extend([fmt, aa.nucleotide])
                    else:
                        segments.append(aa.nucleotide)

                if segments and not isinstance(segments[0], str):
                    segments.insert(0, "")

                worksheet.write_rich_string(row_idx + 1, col_idx, *segments)
            else:
                worksheet.write(row_idx + 1, col_idx, value)

    workbook.close()


def main(
    peptides_path,
    sample_name,
    classI_ic50_score_max,
    classI_ic50_percentile_max,
    classII_ic50_score_max,
    classII_ic50_percentile_max,
    problematic_position,
    output_file,
    output_path,
):

    peptides_df = pd.read_excel(peptides_path)
    peptides_df["RESTRICTING HLA ALLELE"] = ""
    peptides_df["Stylized Sequence"] = pd.Series(
        [[] for _ in range(len(peptides_df))], dtype=object
    )

    for index, row in peptides_df.iterrows():
        restricting_alleles = ""
        if (
            float(row["Class I IC50 MT"]) < classI_ic50_score_max
            or float(row["Class I %ile MT"]) < classI_ic50_percentile_max
        ):
            restricting_alleles = row["Class I Allele"]
        if (
            float(row["Class II IC50 MT"]) < classII_ic50_score_max
            or float(row["Class II %ile MT"]) < classII_ic50_percentile_max
        ):
            if restricting_alleles:
                restricting_alleles += "/" + row["Class II Allele"]
            else:
                restricting_alleles = row["Class II Allele"]

        peptides_df.at[index, "RESTRICTING HLA ALLELE"] = restricting_alleles

        sequence = row[
            "CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES"
        ]
        classI_peptide = row["Best Peptide Class I"]
        classI_ic50 = row["Class I IC50 MT"]
        classI_percentile = row["Class I %ile MT"]
        classI_transcript = row["Class I Best Transcript"]
        classII_peptide = row["Best Peptide Class II"]
        classII_ic50 = row["Class II IC50 MT"]
        classII_percentile = row["Class II %ile MT"]
        classII_transcript = row["Class II Best Transcript"]
        mutant_peptide_pos = str(row["Pos"])

        peptide_sequence = annotate_every_nucleotide(
            sequence,
            classI_peptide,
            classII_peptide,
            classI_ic50,
            classI_percentile,
            classII_ic50,
            classII_percentile,
            classI_transcript,
            classII_transcript,
            classI_ic50_score_max,
            classI_ic50_percentile_max,
            classII_ic50_score_max,
            classII_ic50_percentile_max,
            problematic_position,
        )

        set_underline(peptide_sequence, mutant_peptide_pos, row["full ID"])
        peptides_df.at[index, "Stylized Sequence"] = peptide_sequence

    generate_formatted_excel(peptides_df, output_path, output_file, sample_name)


if __name__ == "__main__":
    main()
