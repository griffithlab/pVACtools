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
    return files[0] if files else None


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
    output_path = os.path.join(output_dir, "mhc_class_i" if class_type == "1" else "mhc_class_ii")

    yml1_path = find_file(results_folder1, os.path.join(folder1_prefix, "log"), "inputs.yml")
    yml2_path = find_file(results_folder2, os.path.join(folder2_prefix, "log"), "inputs.yml")
    if not yml1_path:
        file_name = "inputs_class_I.yml" if class_type == 1 else "inputs_class_II.yml"
        yml1_path = find_file(results_folder1, os.path.dirname(folder1_prefix), file_name)
    if not yml2_path:
        file_name = "inputs_class_I.yml" if class_type == 1 else "inputs_class_II.yml"
        yml2_path = find_file(results_folder2, os.path.dirname(folder2_prefix), file_name)
    if yml1_path and yml2_path:
        logging.info("Running the input YML comparison tool...")
        run_compare_yml(yml1_path, yml2_path, output_path, class_type)
        logging.info("\u2713 Comparison completed successfully.")
    else:
        if not yml1_path and not yml2_path:
            logging.error(
                "ERROR: Could not locate the input YML file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not yml1_path:
            logging.error(
                "ERROR: Could not locate the input YML file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not yml2_path:
            logging.error(
                "ERROR: Could not locate the input YML file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        logging.info("\u2716 Comparison skipped.")

    json1_path = find_file(
        results_folder1, folder1_prefix, "*all_epitopes.aggregated.metrics.json"
    )
    json2_path = find_file(
        results_folder2, folder2_prefix, "*all_epitopes.aggregated.metrics.json"
    )
    if json1_path and json2_path:
        logging.info("\nRunning the metrics JSON comparison tool...")
        run_compare_json(json1_path, json2_path, output_path, class_type)
        logging.info("\u2713 Comparison completed successfully.")
    else:
        if not json1_path and not json2_path:
            logging.error(
                "ERROR: Could not locate the metrics JSON file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not json1_path:
            logging.error(
                "ERROR: Could not locate the metrics JSON file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not json2_path:
            logging.error(
                "ERROR: Could not locate the metrics JSON file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        logging.info("\u2716 Comparison skipped.")

    agg_tsv1_path = find_file(
        results_folder1, folder1_prefix, "*all_epitopes.aggregated.tsv"
    )
    agg_tsv2_path = find_file(
        results_folder2, folder2_prefix, "*all_epitopes.aggregated.tsv"
    )
    if agg_tsv1_path and agg_tsv2_path:
        logging.info("\nRunning the aggregated TSV comparison tool...")
        run_compare_aggregated_tsv(
            agg_tsv1_path, agg_tsv2_path, aggregated_columns, output_path, class_type
        )
        logging.info("\u2713 Comparison completed successfully.")
    else:
        if not agg_tsv1_path and not agg_tsv2_path:
            logging.error(
                "ERROR: Could not locate the aggregated TSV file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not agg_tsv1_path:
            logging.error(
                "ERROR: Could not locate the aggregated TSV file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not agg_tsv2_path:
            logging.error(
                "ERROR: Could not locate the aggregated TSV file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        logging.info("\u2716 Comparison skipped.")

    unagg_tsv1_path = find_file(
        results_folder1, folder1_prefix, "*all_epitopes.tsv"
    )
    unagg_tsv2_path = find_file(
        results_folder2, folder2_prefix, "*all_epitopes.tsv"
    )
    if unagg_tsv1_path and unagg_tsv2_path:
        logging.info("\nRunning the unaggregated TSV comparison tool...")
        run_compare_unaggregated_tsv(
            unagg_tsv1_path,
            unagg_tsv2_path,
            unaggregated_columns,
            output_path,
            class_type,
        )
        logging.info("\u2713 Comparison completed successfully.")
    else:
        if not unagg_tsv1_path and not unagg_tsv2_path:
            logging.error(
                "ERROR: Could not locate the unaggregated TSV file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not unagg_tsv1_path:
            logging.error(
                "ERROR: Could not locate the unaggregated TSV file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not unagg_tsv2_path:
            logging.error(
                "ERROR: Could not locate the unaggregated TSV file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        logging.info("\u2716 Comparison skipped.")

    refmatch_tsv1_path = find_file(
        results_folder1, folder1_prefix, "*.reference_matches"
    )
    refmatch_tsv2_path = find_file(
        results_folder2, folder2_prefix, "*.reference_matches"
    )
    if refmatch_tsv1_path and refmatch_tsv2_path:
        logging.info("\nRunning the reference match TSV comparison tool...")
        run_compare_reference_matches_tsv(
            refmatch_tsv1_path,
            refmatch_tsv2_path,
            reference_match_columns,
            output_path,
            class_type,
        )
        logging.info("\u2713 Comparison completed successfully.")
    else:
        if not refmatch_tsv1_path and not refmatch_tsv2_path:
            logging.error(
                "ERROR: Could not locate the reference match TSV file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not refmatch_tsv1_path:
            logging.error(
                "ERROR: Could not locate the reference match TSV file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif not refmatch_tsv2_path:
            logging.error(
                "ERROR: Could not locate the reference match TSV file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        logging.info("\u2716 Comparison skipped.")
    logging.info("\n" + "\u2500" * 55)
    logging.info(
        "Successfully generated MHC Class %s comparison report.",
        "I" if class_type == "1" else "II",
    )
    logging.info("\u2500" * 55)
