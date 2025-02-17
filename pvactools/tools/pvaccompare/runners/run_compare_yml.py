import logging
from ..run_utils import *
from ..comparisons import CompareYML


def main(input_file1, input_file2, output_path, class_type):
    """
    Purpose:    Control function for the inputs.yml file comparison
    Modifies:   Nothing
    Returns:    None
    """
    logging.basicConfig(level=logging.INFO)
    comparer = CompareYML(input_file1, input_file2)

    try:
        differences = comparer.convert_diff_to_dict()
        if not differences:
            logging.info("The YAML input files are identical.")
    except Exception as e:
        logging.error(f"Error occurred while generating input comparison report: {e}")

    export_to_json(
        comparer.input_file1,
        comparer.input_file2,
        differences,
        "yml_input_data.json",
        output_path,
        class_type,
    )


if __name__ == "__main__":
    import sys

    main(sys.argv[1:])
