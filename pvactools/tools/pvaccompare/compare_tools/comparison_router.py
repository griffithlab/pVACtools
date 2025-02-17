import glob
import os
import logging
from ..runners import *


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
    if os.path.exists(os.path.join(results_folder, "MHC_Class_I")) or os.path.exists(
        os.path.join(results_folder, "MHC_Class_II")
    ):
        return "MHC_Class_I" if class_type == "1" else "MHC_Class_II"
    elif os.path.exists(
        os.path.join(results_folder, "pVACseq/mhc_i")
    ) or os.path.exists(os.path.join(results_folder, "pVACseq/mhc_ii")):
        return "pVACseq/mhc_i" if class_type == "1" else "pVACseq/mhc_ii"
    else:
        raise FileNotFoundError(
            f"Could not locate result files for folder: {results_folder}"
        )


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
    output_path = (
        f'{output_dir}/{"mhc_class_i" if class_type == "1" else "mhc_class_ii"}'
    )

    yml1_path = find_file(results_folder1, folder1_prefix + "/log", "inputs.yml")
    yml2_path = find_file(results_folder2, folder2_prefix + "/log", "inputs.yml")
    if yml1_path and yml2_path:
        logging.info("Running the input YML comparison tool...")
        run_compare_yml(yml1_path, yml2_path, output_path, class_type)
        logging.info("\u2713 Comparison completed successfully.")
    else:
        if yml1_path:
            logging.error(
                "ERROR: Could not locate the input YML file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif yml2_path:
            logging.error(
                "ERROR: Could not locate the input YML file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        else:
            logging.error(
                "ERROR: Could not locate the input YML file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )

        logging.info("\u2716 Comparison skipped.")

    json1_path = find_file(
        results_folder1, folder1_prefix + "/", "*all_epitopes.aggregated.metrics.json"
    )
    json2_path = find_file(
        results_folder2, folder2_prefix + "/", "*all_epitopes.aggregated.metrics.json"
    )
    if json1_path and json2_path:
        logging.info("\nRunning the metrics JSON comparison tool...")
        run_compare_json(json1_path, json2_path, output_path, class_type)
        logging.info("\u2713 Comparison completed successfully.")
    else:
        if json1_path:
            logging.error(
                "ERROR: Could not locate the metrics JSON file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif json2_path:
            logging.error(
                "ERROR: Could not locate the metrics JSON file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        else:
            logging.error(
                "ERROR: Could not locate the metrics JSON file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        logging.info("\u2716 Comparison skipped.")

    agg_tsv1_path = find_file(
        results_folder1, folder1_prefix + "/", "*all_epitopes.aggregated.tsv"
    )
    agg_tsv2_path = find_file(
        results_folder2, folder2_prefix + "/", "*all_epitopes.aggregated.tsv"
    )
    if agg_tsv1_path and agg_tsv2_path:
        logging.info("\nRunning the aggregated TSV comparison tool...")
        run_compare_aggregated_tsv(
            agg_tsv1_path, agg_tsv2_path, aggregated_columns, output_path, class_type
        )
        logging.info("\u2713 Comparison completed successfully.")
    else:
        if agg_tsv1_path:
            logging.error(
                "ERROR: Could not locate the aggregated TSV file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif agg_tsv2_path:
            logging.error(
                "ERROR: Could not locate the aggregated TSV file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        else:
            logging.error(
                "ERROR: Could not locate the aggregated TSV file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        logging.info("\u2716 Comparison skipped.")

    unagg_tsv1_path = find_file(
        results_folder1, folder1_prefix + "/", "*all_epitopes.tsv"
    )
    unagg_tsv2_path = find_file(
        results_folder2, folder2_prefix + "/", "*all_epitopes.tsv"
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
        if unagg_tsv1_path:
            logging.error(
                "ERROR: Could not locate the unaggregated TSV file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif unagg_tsv2_path:
            logging.error(
                "ERROR: Could not locate the unaggregated TSV file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        else:
            logging.error(
                "ERROR: Could not locate the unaggregated TSV file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        logging.info("\u2716 Comparison skipped.")

    refmatch_tsv1_path = find_file(
        results_folder1, folder1_prefix + "/", "*.reference_matches"
    )
    refmatch_tsv2_path = find_file(
        results_folder2, folder2_prefix + "/", "*.reference_matches"
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
        if refmatch_tsv1_path:
            logging.error(
                "ERROR: Could not locate the reference match TSV file in results folder 2 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        elif refmatch_tsv2_path:
            logging.error(
                "ERROR: Could not locate the reference match TSV file in results folder 1 for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        else:
            logging.error(
                "ERROR: Could not locate the reference match TSV file in either results folder for MHC Class %s.",
                "I" if class_type == "1" else "II",
            )
        logging.info("\u2716 Comparison skipped.")
    logging.info("\n" + "\u2500" * 55)
    logging.info(
        "Successfully generated MHC Class %s comparison report.",
        "I" if class_type == "1" else "II",
    )
    logging.info("\u2500" * 55)
