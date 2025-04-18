import pandas as pd
import numpy as np
import re
import logging
import json


def add_line_numbers(df1, df2):
    """
    Purpose:    Add line number information to each row of the given DataFrames.
    Modifies:   df1 and df2
    Returns:    None
    """
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


def find_dropped_cols(df1, df2, original_columns, id_columns):
    """
    Purpose:    Outputs the dropped comparison columns to the terminal and creates an array of notes for the html report
    Modifies:   Nothing
    Returns:    List of strings, run_notes
    """
    run_notes = []
    for col in original_columns:
        prefix = "ID " if col in id_columns else ""
        if col not in df1.columns and col not in df2.columns:
            message = (
                f"\u2022 {prefix}Column dropped: '{col}' is not present in either file"
            )
        elif col not in df1.columns:
            message = (
                f"\u2022 {prefix}Column dropped: '{col}' is only present in file 2"
            )
        elif col not in df2.columns:
            message = (
                f"\u2022 {prefix}Column dropped: '{col}' is only present in file 1"
            )
        else:
            continue
        logging.info(message)
        run_notes.append(message.replace("\u2022 ", ""))
    return run_notes


def get_common_entries(df1, df2):
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


def get_unique_entries(df1, df2, common_entries, contains_id=True):
    """
    Purpose:    Find, store, and sort unique entries to each dataframe
    Modifies:   Nothing
    Returns:    Two sets containing IDs unique to the corresponding dataframes
    """
    return sort_unique_entries(
        set(df1["ID"]).difference(common_entries),
        set(df2["ID"]).difference(common_entries),
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
        logging.error(f"ERROR: While splitting replaced ID: {id_str}, {e}")
        return "", ""


def sort_unique_entries(unique_entries_file1, unique_entries_file2, contains_id):
    """
    Purpose:    Sort unique entries by ID or Gene-AA combination depending on input
    Modifies:   Nothing
    Returns:    Two sorted lists of unique entries from file 1 and file 2
    """
    if unique_entries_file1:
        if contains_id:
            unique_entries_file1 = sorted(unique_entries_file1, key=extract_id_parts)
        else:
            unique_entries_file1 = sorted(unique_entries_file1, key=split_replaced_id)
    if unique_entries_file2:
        if contains_id:
            unique_entries_file2 = sorted(unique_entries_file2, key=extract_id_parts)
        else:
            unique_entries_file2 = sorted(unique_entries_file2, key=split_replaced_id)
    return unique_entries_file1, unique_entries_file2


def get_file_differences(
    df1,
    df2,
    id_columns,
    columns_to_compare,
    contains_id=True,
    tolerance=0.1,
):
    """
    Purpose:    Find and store differences found between the two dataframes
    Modifies:   Nothing
    Returns:    Dictionary of differences and a dictionary of unique entries
    """
    id_cols = ["ID", "line"] + id_columns
    df1_selected = df1[id_cols + columns_to_compare]
    df2_selected = df2[id_cols + columns_to_compare]

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

        keep_cols = ["ID", col_file1, col_file2, "line_file1", "line_file2"]
        for colname in id_columns:
            if f"{colname}_file1" in merged_df.columns:
                keep_cols.append(f"{colname}_file1")

        diff = merged_df[mask][keep_cols]
        diff = diff.rename(columns={f"{col}_file1": col for col in id_columns})

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
    """
    Purpose:    Recursively replace NaN float values with None in nested lists or dictionaries
    Modifies:   Nothing
    Returns:    The same object with all NaNs replaced by None
    """
    if isinstance(obj, float) and pd.isna(obj):
        return None
    elif isinstance(obj, dict):
        return {k: replace_nan_with_none(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [replace_nan_with_none(v) for v in obj]
    return obj


def merge_hla_field(split_id_parts):
    """
    Purpose: Combines 'HLA' and the following part into a single element if found
    Modifies: Nothing
    Returns: List of properly split ID parts
    """
    new_parts = []
    skip = False
    for i, part in enumerate(split_id_parts):
        if skip:
            skip = False
            continue
        if part == "HLA" and i + 1 < len(split_id_parts):
            new_parts.append(f"HLA-{split_id_parts[i+1]}")
            skip = True
        else:
            new_parts.append(part)
    return new_parts


def fill_entry_dict(split_id, id_columns, entry_num, key):
    """
    Purpose:    Build a dictionary that maps split ID parts to their corresponding ID columns
    Modifies:   Nothing
    Returns:    A dictionary with ID column values and entry number, or None if mismatch occurs
    """
    if len(split_id) != len(id_columns):
        logging.error(
            f"ERROR: Mismatch between number of ID columns and split ID for key: {key}"
        )
        return None

    entry_dict = {"Entry #": entry_num}
    for i, col in enumerate(id_columns):
        entry_dict[col] = split_id[i]
    return entry_dict


def preprocess_unique_entries(entries, id_columns):
    """
    Purpose:    Convert unique entry sets or dictionaries into a structured list of dictionaries
    Modifies:   Nothing
    Returns:    A list of structured entry dictionaries with ID fields and optional hit counts
    """
    transformed_entries = []
    entry_num = 1

    if isinstance(entries, dict):
        for key, value in entries.items():
            split_id = merge_hla_field(key.split("-"))
            entry_dict = fill_entry_dict(split_id, id_columns, entry_num, key)
            if entry_dict:
                entry_dict["Number of Hits"] = value
                transformed_entries.append(entry_dict)
                entry_num += 1
    elif isinstance(entries, list):
        for entry_id in entries:
            split_id = merge_hla_field(entry_id.split("-"))
            entry_dict = fill_entry_dict(split_id, id_columns, entry_num, entry_id)
            if entry_dict:
                transformed_entries.append(entry_dict)
                entry_num += 1
    else:
        logging.error(
            "ERROR: Could not process unique entries. Must be a dictionary or a list."
        )
    return transformed_entries


def preprocess_differences(differences, id_columns, chunk_size=1000):
    """
    Purpose:    Chunk and restructure difference records for display/export, replacing NaNs with None
    Modifies:   Nothing
    Returns:    A dictionary of chunked, structured difference entries with metadata
    """
    differences = replace_nan_with_none(differences)
    transformed_differences = {}

    for section, entries in differences.items():
        entry_num = 1
        chunked_entries = [
            entries[i : i + chunk_size] for i in range(0, len(entries), chunk_size)
        ]

        section_data = {"num_sections": len(chunked_entries)}
        section_content = {}

        for idx, chunk in enumerate(chunked_entries):
            section_content[f"section{idx + 1}"] = []

            for i, entry in enumerate(chunk):
                entry_dict = {
                    "Entry #": entry_num + i,
                }

                for col in id_columns:
                    entry_dict[col] = entry.get(col)

                entry_dict.update(
                    {
                        "File 1 Value": entry.get(f"{section}_file1"),
                        "File 2 Value": entry.get(f"{section}_file2"),
                        "File 1 Line": entry.get("line_file1"),
                        "File 2 Line": entry.get("line_file2"),
                    }
                )

                section_content[f"section{idx + 1}"].append(entry_dict)
            entry_num += len(chunk)

        section_data.update(section_content)
        transformed_differences[section] = section_data

    return transformed_differences


def export_to_json(
    input_file1,
    input_file2,
    differences,
    filename,
    output_path,
    class_type,
    id_columns=[],
    run_notes=[],
    common_entries=[],
    unique_entries_file1=[],
    unique_entries_file2=[],
    hits_file1={},
    hits_file2={},
    duplicate_ids=False,
):
    """
    Purpose: Exports comparison results and statistics between two input files into a JSON file
    Modifies: Writes a structured JSON file to disk at the specified output path
    Returns: None
    """
    file_path = f"{output_path}/{filename}"

    if filename != "yml_input_data.json" and filename != "json_input_data.json":
        summary_data = {
            "Notes": run_notes,
            "Entries": {
                "Total number of entries": get_total_number_entries(
                    common_entries, unique_entries_file1, unique_entries_file2
                ),
                "Number of common entries": len(common_entries),
                "Number of entries unique to file 1": len(unique_entries_file1),
                "Number of entries unique to file 2": len(unique_entries_file2),
            },
            "Section Differences": {},
        }
        num_col_differences = get_number_column_differences(differences)
        for col, _ in differences.items():
            if col != "ID":
                summary_data["Section Differences"][
                    f"Number of differences in {col}"
                ] = num_col_differences[col]
        differences = preprocess_differences(differences, id_columns)
    else:
        summary_data = {}

    entry_data = {
        "Entries Unique to File 1" if not duplicate_ids else "Hits in File 1": (
            preprocess_unique_entries(list(unique_entries_file1), id_columns)
            if not duplicate_ids
            else preprocess_unique_entries(hits_file1, id_columns)
        ),
        "Entries Unique to File 2" if not duplicate_ids else "Hits in File 2": (
            preprocess_unique_entries(list(unique_entries_file2), id_columns)
            if not duplicate_ids
            else preprocess_unique_entries(hits_file2, id_columns)
        ),
    }

    data = {
        "mhc_class": class_type,
        "input_file1": input_file1,
        "input_file2": input_file2,
        "summary": summary_data,
        "differences": differences,
        "entries": entry_data,
    }

    with open(file_path, "w") as f:
        json.dump(data, f, indent=4)


def get_total_number_entries(
    common_entries, unique_entries_file1, unique_entries_file2
):
    """
    Purpose:    Get the total number of entries between the two files
    Modifies:   Nothing
    Returns:    Integer of the total number of entries
    """
    return len(common_entries) + len(unique_entries_file1) + len(unique_entries_file2)


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


def check_id_columns(df1, df2, id_cols):
    """
    Purpose:    Makes sure the 'ID' columns are present in both dataframes
    Modifies:   None
    Returns:    Updated 'ID' column list
    """
    updated_id_cols = []
    for col in id_cols:
        if col in df1.columns and col in df2.columns:
            updated_id_cols.append(col)
    return updated_id_cols


def create_id_column(df1, df2, cols):
    """
    Purpose:    Combines multiple columns into a singular unique ID column in both dataframes
    Modifies:   df1 and df2
    Returns:    Modified dataframes
    """
    df1["ID"] = df1[cols].apply(lambda x: "-".join(map(str, x)), axis=1)
    df2["ID"] = df2[cols].apply(lambda x: "-".join(map(str, x)), axis=1)
