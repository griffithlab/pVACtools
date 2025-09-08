import glob
import os
import logging
from pvactools.tools.pvaccompare.runners import *


def find_file(results_folder, subfolder, pattern):
    """
    Purpose:    Attempts to locate the files needed for each comparison
    Modifies:   Nothing
    Returns:    A string of the file path
    """
    search_path = os.path.join(results_folder, subfolder, pattern)
    files = glob.glob(search_path, recursive=True)
    return files


def get_prefix(class_type, results_folder):
    results_folder = os.path.abspath(results_folder)
    matches = []

    # Walk through all subdirectories recursively
    for root, dirs, _ in os.walk(results_folder):
        for d in dirs:
            d_lower = d.lower()
            if class_type == "1":
                if "mhc" in d_lower and "i" in d_lower and "ii" not in d_lower:
                    matches.append(os.path.join(root, d))
            elif class_type == "2":
                if "mhc" in d_lower and "ii" in d_lower:
                    matches.append(os.path.join(root, d))

    if not matches:
        raise FileNotFoundError(
            f"Could not locate result files for folder: {results_folder} and class {class_type}"
        )
    if len(matches) > 1:
        raise RuntimeError(
            f"Multiple matching directories found for class {class_type}: {matches}"
        )

    return matches[0]


def run_comparison(
    class_type, file1_path, file2_path, description="", run_func=None, run_args=()
):
    mhc_class = "I" if class_type == "1" else "II"

    if file1_path and file2_path and len(file1_path) == 1 and len(file2_path) == 1:
        logging.info(f"Running the {description} comparison tool...")
        run_func(file1_path[0], file2_path[0], *run_args)
        logging.info("\u2713 Comparison completed successfully.\n")
        return

    error_conditions = [
        (
            not file1_path,
            "Could not locate the %s file in results folder 1 for MHC Class %s.",
            None,
        ),
        (
            not file2_path,
            "Could not locate the %s file in results folder 2 for MHC Class %s.",
            None,
        ),
        (
            file1_path and len(file1_path) > 1,
            "Located multiple %s files in results folder 1 for MHC Class %s:",
            file1_path,
        ),
        (
            file2_path and len(file2_path) > 1,
            "Located multiple %s files in results folder 2 for MHC Class %s:",
            file2_path,
        ),
    ]

    for condition, message, extra in error_conditions:
        if condition:
            logging.error("ERROR: " + message, description, mhc_class)
            if extra:  # extra contains the file paths if multiple were found
                for f in extra:
                    logging.error("     %s", f)

    logging.info("\u2716 Comparison skipped.\n")


def main(
    class_type,
    results_folder1,
    results_folder2,
    output_dir,
    aggregated_columns,
    unaggregated_columns,
    reference_match_columns,
):
    """
    Purpose:    Runs all of the different comparisons
    Modifies:   Nothing
    Returns:    None
    """
    folder1_prefix = get_prefix(class_type, results_folder1)
    folder2_prefix = get_prefix(class_type, results_folder2)
    output_path = os.path.join(
        output_dir, "mhc_class_i" if class_type == "1" else "mhc_class_ii"
    )

    # Input YML comparison
    yml1_path = find_file(
        results_folder1, os.path.join(folder1_prefix, "log"), "inputs.yml"
    )
    yml2_path = find_file(
        results_folder2, os.path.join(folder2_prefix, "log"), "inputs.yml"
    )
    if not yml1_path:
        file_name = "inputs_class_I.yml" if class_type == 1 else "inputs_class_II.yml"
        yml1_path = find_file(
            results_folder1, os.path.dirname(folder1_prefix), file_name
        )
    if not yml2_path:
        file_name = "inputs_class_I.yml" if class_type == 1 else "inputs_class_II.yml"
        yml2_path = find_file(
            results_folder2, os.path.dirname(folder2_prefix), file_name
        )
    run_comparison(
        class_type=class_type,
        file1_path=yml1_path,
        file2_path=yml2_path,
        description="input YML",
        run_func=run_compare_yml,
        run_args=(output_path, class_type),
    )

    # Metrics JSON comparison
    json1_path = find_file(
        results_folder1, folder1_prefix, "*all_epitopes.aggregated.metrics.json"
    )
    json2_path = find_file(
        results_folder2, folder2_prefix, "*all_epitopes.aggregated.metrics.json"
    )
    run_comparison(
        class_type=class_type,
        file1_path=json1_path,
        file2_path=json2_path,
        description="metrics JSON",
        run_func=run_compare_json,
        run_args=(output_path, class_type),
    )

    # Aggregated TSV comparison
    agg_tsv1_path = find_file(
        results_folder1, folder1_prefix, "*all_epitopes.aggregated.tsv"
    )
    agg_tsv2_path = find_file(
        results_folder2, folder2_prefix, "*all_epitopes.aggregated.tsv"
    )
    run_comparison(
        class_type=class_type,
        file1_path=agg_tsv1_path,
        file2_path=agg_tsv2_path,
        description="aggregated TSV",
        run_func=run_compare_aggregated_tsv,
        run_args=(aggregated_columns, output_path, class_type),
    )

    # Unaggregated TSV comparison
    unagg_tsv1_path = find_file(results_folder1, folder1_prefix, "*all_epitopes.tsv")
    unagg_tsv2_path = find_file(results_folder2, folder2_prefix, "*all_epitopes.tsv")
    run_comparison(
        class_type=class_type,
        file1_path=unagg_tsv1_path,
        file2_path=unagg_tsv2_path,
        description="unaggregated TSV",
        run_func=run_compare_unaggregated_tsv,
        run_args=(unaggregated_columns, output_path, class_type),
    )

    # Reference matches TSV comparison
    refmatch_tsv1_path = find_file(
        results_folder1, folder1_prefix, "*.reference_matches"
    )
    refmatch_tsv2_path = find_file(
        results_folder2, folder2_prefix, "*.reference_matches"
    )
    run_comparison(
        class_type=class_type,
        file1_path=refmatch_tsv1_path,
        file2_path=refmatch_tsv2_path,
        description="reference match TSV",
        run_func=run_compare_reference_matches_tsv,
        run_args=(reference_match_columns, output_path, class_type),
    )

    logging.info("\u2500" * 55)
    logging.info(
        "Successfully generated MHC Class %s comparison report.",
        "I" if class_type == "1" else "II",
    )
    logging.info("\u2500" * 55 + "\n")
