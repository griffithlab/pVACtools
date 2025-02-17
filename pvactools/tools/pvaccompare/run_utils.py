import pandas as pd
import numpy as np
import re
import logging
import json


def add_line_numbers(df1, df2):
    df1["line"] = range(2, len(df1) + 2)
    df2["line"] = range(2, len(df2) + 2)


def check_column_formatting(df1, df2):
    """
    Purpose:    Rename columns based on the mappings dictionary to make column names the same
    Modifies:   df1 and df2
    Returns:    None
    """
    column_mappings = {  # Fill in different names/formatting between versions
        "Best Peptide": ["best peptide", "best_peptide"],
        "Best Transcript": ["best transcript", "best_transcript"],
        "Tier": ["tier"],
        "AA Change": ["AA_change"],
        "Num Passing Transcripts": ["Num_Transcript"],
        "Num Passing Peptides": ["Num_Peptides"],
    }

    for col in df1.columns:
        for key, value in column_mappings.items():
            if col == key:
                break
            elif col in value:
                logging.info("\u2022 Renamed '%s' to '%s' in file 1", col, key)
                df1.rename(columns={col: key}, inplace=True)
                break
    for col in df2.columns:
        for key, value in column_mappings.items():
            if col == key:
                break
            elif col in value:
                logging.info("\u2022 Renamed '%s' to '%s' in file 2", col, key)
                df2.rename(columns={col: key}, inplace=True)
                break


def find_dropped_cols(df1, df2, original_columns):
    """
    Purpose:    Outputs the dropped comparison columns to the terminal and creates an array of notes for the html report
    Modifies:   Nothing
    Returns:    List of strings, run_notes
    """
    run_notes = []
    for col in original_columns:
        if col not in df1.columns and col not in df2.columns:
            logging.info(
                "\u2022 Column dropped: '%s' is not present in either file", col
            )
            run_notes.append(f"Column dropped: '{col}' is not present in either file")
        elif col not in df1.columns:
            logging.info("\u2022 Column dropped: '%s' is only present in file 2", col)
            run_notes.append(f"Column dropped: '{col}' is only present in file 2")
        elif col not in df2.columns:
            logging.info("\u2022 Column dropped: '%s' is only present in file 1", col)
            run_notes.append(f"Column dropped: '{col}' is only present in file 1")
    return run_notes


def get_common_variants(df1, df2):
    """
    Purpose:    Find and store IDs shared between the two given dataframes
    Modifies:   Nothing
    Returns:    A set containing IDs that are common between the two dataframes
    """
    return set(df1["ID"]).intersection(set(df2["ID"]))


def load_tsv_files(input_file1, input_file2):
    """
    Purpose:    Load the two input tsv files into dataframes
    Modifies:   Nothing
    Returns:    Two dataframes corresponding to the two input files
    """
    try:
        df1 = pd.read_csv(input_file1, sep="\t", low_memory=False)
        df2 = pd.read_csv(input_file2, sep="\t", low_memory=False)
    except Exception as e:
        raise Exception(f"Error loading files: {e}")
    return df1, df2


def check_columns_to_compare(df1, df2, columns_to_compare):
    """
    Purpose:    Add columns present in both dataframes to columns_to_keep
    Modifies:   Nothing
    Returns:    List of columns present in both dataframes
    """
    columns_to_keep = []
    for col in columns_to_compare:
        if col in df1.columns and col in df2.columns:
            columns_to_keep.append(col)
    return columns_to_keep


def get_unique_variants(df1, df2, common_variants, contains_id=True):
    """
    Purpose:    Find, store, and sort unique variants to each dataframe
    Modifies:   Nothing
    Returns:    Two sets containing IDs unique to the corresponding dataframes
    """
    return sort_unique_variants(
        set(df1["ID"]).difference(common_variants),
        set(df2["ID"]).difference(common_variants),
        contains_id,
    )


def extract_id_parts(id_str):
    """
    Purpose:    Extract parts of the ID to use in sorting
    Modifies:   Nothing
    Returns:    A tuple of the different sections for sorting
    """
    match = re.match(r"chr(\w+)-(\d+)-(\d+)-", id_str)
    if match:
        chr_part = match.group(1)
        if chr_part.isdigit():
            chr_part = int(chr_part)
        else:
            chr_part = float("inf")
        return chr_part, int(match.group(2)), int(match.group(3))
    return None, None, None


def split_replaced_id(id_str):
    """
    Purpose:    Extract parts of the replaced ID (Gene-AA change) to use in sorting
    Modifies:   Nothing
    Returns:    Two strings corresponding to the split sections
    """
    try:
        grp1, rest = id_str.split(" (")
        grp2 = rest.split("-")[0].rstrip(")")
        return grp1, grp2
    except Exception as e:
        logging.error(f"Error splitting replaced ID: {id_str}, {e}")
        return "", ""


def sort_unique_variants(unique_variants_file1, unique_variants_file2, contains_id):
    if unique_variants_file1:
        if contains_id:
            unique_variants_file1 = sorted(unique_variants_file1, key=extract_id_parts)
        else:
            unique_variants_file1 = sorted(unique_variants_file1, key=split_replaced_id)
    if unique_variants_file2:
        if contains_id:
            unique_variants_file2 = sorted(unique_variants_file2, key=extract_id_parts)
        else:
            unique_variants_file2 = sorted(unique_variants_file2, key=split_replaced_id)
    return unique_variants_file1, unique_variants_file2


def get_file_differences(
    df1,
    df2,
    columns_to_compare,
    contains_id=True,
    tolerance=0.1,
):
    """
    Purpose:    Find and store differences found between the two dataframes
    Modifies:   Nothing
    Returns:    Dictionary of differences and a dictionary of unique variants
    """
    df1_selected = df1[["ID", "line"] + columns_to_compare]
    df2_selected = df2[["ID", "line"] + columns_to_compare]

    merged_df = pd.merge(
        df1_selected, df2_selected, on="ID", suffixes=("_file1", "_file2")
    )

    differences = {}
    for col in columns_to_compare:
        col_file1 = f"{col}_file1"
        col_file2 = f"{col}_file2"

        # Determine if columns are numeric
        is_numeric_col1 = np.issubdtype(merged_df[col_file1].dtype, np.number)
        is_numeric_col2 = np.issubdtype(merged_df[col_file2].dtype, np.number)

        if is_numeric_col1 and is_numeric_col2:
            # Convert to numeric values to ensure uniformity and handle missing/invalid values
            merged_df[col_file1] = pd.to_numeric(merged_df[col_file1], errors="coerce")
            merged_df[col_file2] = pd.to_numeric(merged_df[col_file2], errors="coerce")

            # Mask for numeric differences greater than tolerance
            tolerance_mask = (
                np.abs(merged_df[col_file1] - merged_df[col_file2]) > tolerance
            )

            # Mask for rows where one value is NaN and the other is not
            nan_mask = (merged_df[col_file1].isna() & ~merged_df[col_file2].isna()) | (
                ~merged_df[col_file1].isna() & merged_df[col_file2].isna()
            )

            # Final mask includes rows with significant numeric differences or NaN-regular number comparisons
            mask = tolerance_mask | nan_mask
        else:
            mask = (merged_df[col_file1] != merged_df[col_file2]) & ~(
                merged_df[col_file1].isna() & merged_df[col_file2].isna()
            )

        diff = merged_df[mask][["ID", col_file1, col_file2, "line_file1", "line_file2"]]
        if not diff.empty:
            differences[col] = diff.to_dict("records")

    for col in differences:
        if contains_id:
            differences[col] = sorted(
                differences[col], key=lambda x: extract_id_parts(x["ID"])
            )
        else:
            differences[col] = sorted(
                differences[col], key=lambda x: split_replaced_id(x["ID"])
            )

    return differences


def replace_nan_with_none(obj):
    if isinstance(obj, float) and pd.isna(obj):
        return None
    elif isinstance(obj, dict):
        return {k: replace_nan_with_none(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [replace_nan_with_none(v) for v in obj]
    return obj


def preprocess_differences(differences, chunk_size=1000):
    differences = replace_nan_with_none(differences)
    transformed_differences = {}

    for section, entries in differences.items():
        chunked_entries = [
            entries[i : i + chunk_size] for i in range(0, len(entries), chunk_size)
        ]

        section_data = {"num_sections": len(chunked_entries)}

        section_data.update(
            {
                f"section{idx + 1}": [
                    {
                        "ID": entry["ID"],
                        "File 1 Value": entry.get(f"{section}_file1"),
                        "File 2 Value": entry.get(f"{section}_file2"),
                        "File 1 Line": entry.get("line_file1"),
                        "File 2 Line": entry.get("line_file2"),
                    }
                    for entry in chunk
                ]
                for idx, chunk in enumerate(chunked_entries)
            }
        )

        transformed_differences[section] = section_data

    return transformed_differences


def export_to_json(
    input_file1,
    input_file2,
    differences,
    filename,
    output_path,
    class_type,
    id_format="",
    run_notes=[],
    common_variants=[],
    unique_variants_file1=[],
    unique_variants_file2=[],
    hits_file1={},
    hits_file2={},
    duplicate_ids=False,
):
    file_path = f"{output_path}/{filename}"

    if filename != "yml_input_data.json" and filename != "json_input_data.json":
        summary_data = {
            "Notes": run_notes,
            "Variants": {
                "Total number of variants": get_total_number_variants(
                    common_variants, unique_variants_file1, unique_variants_file2
                ),
                "Number of common variants": len(common_variants),
                "Number of variants unique to file 1": len(unique_variants_file1),
                "Number of variants unique to file 2": len(unique_variants_file2),
            },
            "Section Differences": {},
        }
        num_col_differences = get_number_column_differences(differences)
        for col, _ in differences.items():
            if col != "ID":
                summary_data["Section Differences"][
                    f"Number of differences in {col}"
                ] = num_col_differences[col]
        differences = preprocess_differences(differences)
    else:
        summary_data = {}

    variant_data = {
        "Variants Unique to File 1" if not duplicate_ids else "Hits in File 1": (
            list(unique_variants_file1) if not duplicate_ids else hits_file1
        ),
        "Variants Unique to File 2" if not duplicate_ids else "Hits in File 2": (
            list(unique_variants_file2) if not duplicate_ids else hits_file2
        ),
    }

    data = {
        "mhc_class": class_type,
        "input_file1": input_file1,
        "input_file2": input_file2,
        "id_format": id_format,
        "summary": summary_data,
        "differences": differences,
        "variants": variant_data,
    }

    with open(file_path, "w") as f:
        json.dump(data, f, indent=4)


def get_total_number_variants(
    common_variants, unique_variants_file1, unique_variants_file2
):
    """
    Purpose:    Get the total number of variants between the two files
    Modifies:   Nothing
    Returns:    Integer of the total number of variants
    """
    return (
        len(common_variants) + len(unique_variants_file1) + len(unique_variants_file2)
    )


def get_number_column_differences(differences):
    """
    Purpose:    Get the number of differences for each column
    Modifies:   Nothing
    Returns:    Dictionary of the columns and corresponding differences
    """
    num_col_differences = {}
    for col, differences in differences.items():
        if col != "ID":
            num_col_differences[col] = len(differences)
    return num_col_differences
