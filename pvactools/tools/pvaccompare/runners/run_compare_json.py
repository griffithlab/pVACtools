import logging
from pvactools.tools.pvaccompare.run_utils import *
from pvactools.tools.pvaccompare.comparisons import CompareJSON


def main(input_file1, input_file2, output_path, class_type):
    """
    Purpose:    Control function for the metrics json file comparison
    Modifies:   Nothing
    Returns:    None
    """
    comparer = CompareJSON(input_file1, input_file2)
    comparer.compare_metric_data()

    if not comparer.differences:
        logging.info("The JSON metric inputs are identical.")

    export_to_json(
        comparer.input_file1,
        comparer.input_file2,
        comparer.differences,
        "json_input_data.json",
        output_path,
        class_type,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
