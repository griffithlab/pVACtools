from pvactools.tools.pvaccompare.run_utils import *
from pvactools.tools.pvaccompare.comparisons import CompareReferenceMatchesTSV
import logging


def main(input_file1, input_file2, columns_to_compare, output_path, class_type):
    """
    Purpose:    Control function for the reference matches tsv comparison
    Modifies:   Nothing
    Returns:    None
    """
    id_format = "Chromosome-Start-Stop-Reference-Variant-Transcript-MT_Epitope_Seq-Hit_ID-Match_Start-Match_Stop"
    duplicate_ids = False

    comparer = CompareReferenceMatchesTSV(input_file1, input_file2, columns_to_compare)
    add_line_numbers(comparer.df1, comparer.df2)
    check_column_formatting(comparer.df1, comparer.df2)
    comparer.create_id_column()

    comparer.run_notes = find_dropped_cols(
        comparer.df1, comparer.df2, comparer.columns_to_compare
    )

    comparer.columns_to_compare = check_columns_to_compare(
        comparer.df1, comparer.df2, comparer.columns_to_compare
    )

    common_variants = get_common_variants(comparer.df1, comparer.df2)
    unique_variants_file1, unique_variants_file2 = get_unique_variants(
        comparer.df1, comparer.df2, common_variants
    )

    if comparer.check_duplicate_ids():
        duplicate_ids = True
        differences = {}
    else:
        differences = get_file_differences(
            comparer.df1, comparer.df2, comparer.columns_to_compare
        )

    if (
        not unique_variants_file1
        and not unique_variants_file2
        and not duplicate_ids
        and not differences
    ):
        logging.info("The Reference Matches TSV files are identical.")

    export_to_json(
        comparer.input_file1,
        comparer.input_file2,
        differences,
        "reference_matches_data.json",
        output_path,
        class_type,
        id_format,
        comparer.run_notes,
        common_variants,
        unique_variants_file1,
        unique_variants_file2,
        comparer.hits_file1,
        comparer.hits_file2,
        duplicate_ids,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
