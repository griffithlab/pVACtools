import pandas as pd
import os
import xlsxwriter
from Bio import SeqIO, pairwise2


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


def get_mutant_positions_from_fasta(fasta_path, full_id):
    if full_id.startswith("WT.") or full_id.startswith("MT."):
        full_id = full_id[3:]

    frameshift = ".FS." in full_id
    wt_id = f"WT.{full_id}"
    mt_id = f"MT.{full_id}"
    wt_seq = None
    mt_seq = None

    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id == wt_id:
            wt_seq = str(record.seq)
        elif record.id == mt_id:
            mt_seq = str(record.seq)

    if wt_seq is None or mt_seq is None:
        print(f"Warning: Missing WT or MT for {full_id}")
        return set()

    # Perform global alignment
    alignments = pairwise2.align.globalms(wt_seq, mt_seq, 2, -1, -5, -0.5)
    wt_aligned, mt_aligned = alignments[0].seqA, alignments[0].seqB

    mutant_positions = set()
    mt_index = 0

    for i in range(len(wt_aligned)):
        wt_residue = wt_aligned[i]
        mt_residue = mt_aligned[i]

        if mt_residue == "-":
            if not frameshift:  # For inframe deletions only, underline residues adjacent to deletions
                if mt_index > 0:
                    mutant_positions.add(mt_index - 1)
                if mt_index < len(mt_seq):
                    mutant_positions.add(mt_index)
        else:
            if wt_residue != mt_residue:
                mutant_positions.add(mt_index)
            mt_index += 1

    return mutant_positions


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


def set_underline(peptide_sequence, mutant_positions):
    for pos in mutant_positions:
        peptide_sequence[pos].underline = True


def generate_formatted_excel(peptides_df, output_path, output_file_prefix, sample_name):
    file_name = f"{output_file_prefix}_{sample_name}.Colored_Peptides.xlsx"
    file_path = os.path.join(output_path, file_name)
    workbook = xlsxwriter.Workbook(file_path)
    worksheet = workbook.add_worksheet()
    format_cache = {}

    def get_format(aa):
        fmt_key = (
            aa.bold,
            aa.color,
            aa.underline,
            aa.large,
        )
        if fmt_key in format_cache:
            return format_cache[fmt_key]

        fmt_props = {}
        if aa.bold:
            fmt_props["bold"] = True
        if aa.color:
            fmt_props["font_color"] = "red"
        if aa.underline:
            fmt_props["underline"] = True
        if aa.large:
            fmt_props["font_size"] = 14

        if fmt_props:
            cell_format = workbook.add_format(fmt_props)
            format_cache[fmt_key] = cell_format
            return cell_format
        else:
            return None

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
                    fmt = get_format(aa)
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
    fasta_path,
    peptides_path,
    sample_name,
    classI_ic50_score_max,
    classI_ic50_percentile_max,
    classII_ic50_score_max,
    classII_ic50_percentile_max,
    problematic_position,
    output_file_prefix,
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

        mutant_positions = get_mutant_positions_from_fasta(fasta_path, row["full ID"])
        set_underline(peptide_sequence, mutant_positions)
        peptides_df.at[index, "Stylized Sequence"] = peptide_sequence

    generate_formatted_excel(peptides_df, output_path, output_file_prefix, sample_name)


if __name__ == "__main__":
    main()
