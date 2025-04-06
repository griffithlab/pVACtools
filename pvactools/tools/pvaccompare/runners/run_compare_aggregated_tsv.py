from pvactools.tools.pvaccompare.comparisons import CompareAggregatedTSV
from pvactools.tools.pvaccompare.run_utils import *
import logging


def main(input_file1, input_file2, columns_to_compare, output_path, class_type):
    """
    Purpose:    Control function for the aggregated tsv file comparison
    Modifies:   Nothing
    Returns:    None
    """
    comparer = CompareAggregatedTSV(input_file1, input_file2, columns_to_compare)
    add_line_numbers(comparer.df1, comparer.df2)
    check_column_formatting(comparer.df1, comparer.df2)
    comparer.check_id()

    id_format = (
        "Chromosome-Start-Stop-Reference-Variant"
        if comparer.contains_id
        else "Gene (AA_Change)"
    )

    run_notes = find_dropped_cols(
        comparer.df1, comparer.df2, comparer.columns_to_compare
    )

    if not comparer.contains_id:
        run_notes.append("Replaced ID with Gene and AA Change")

    comparer.columns_to_compare = check_columns_to_compare(
        comparer.df1, comparer.df2, comparer.columns_to_compare
    )

    common_entries = get_common_entries(comparer.df1, comparer.df2)
    unique_entries_file1, unique_entries_file2 = get_unique_entries(
        comparer.df1, comparer.df2, common_entries, comparer.contains_id
    )

    differences = get_file_differences(
        comparer.df1,
        comparer.df2,
        comparer.columns_to_compare,
        comparer.contains_id,
    )

    if not unique_entries_file1 and not unique_entries_file2 and not differences:
        logging.info("The Aggregated TSV files are identical.")

    export_to_json(
        comparer.input_file1,
        comparer.input_file2,
        differences,
        "aggregated_data.json",
        output_path,
        class_type,
        id_format,
        run_notes,
        common_entries,
        unique_entries_file1,
        unique_entries_file2,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
