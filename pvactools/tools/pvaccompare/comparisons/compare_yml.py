import yaml
from deepdiff import DeepDiff
import re


class CompareYML:
    def __init__(self, input_file1, input_file2):
        self.input_file1 = input_file1
        self.input_file2 = input_file2
        self.data1, self.data2 = self.load_files()
        self.differences = DeepDiff(self.data1, self.data2, ignore_order=True)
        self.output_mappings = {
            "dictionary_item_added": "Fields Unique to File 2",
            "dictionary_item_removed": "Fields Unique to File 1",
            "values_changed": "Values Changed",
            "iterable_item_added": "Values Added in File 2",
            "iterable_item_removed": "Values Removed in File 2",
        }

    def load_files(self):
        """
        Purpose:    Load the two input yml files into dictionaries
        Modifies:   Nothing
        Returns:    Two dictionaries corresponding to the two input files
        """
        with open(self.input_file1, "r") as f1, open(self.input_file2, "r") as f2:
            data1 = yaml.safe_load(f1)
            data2 = yaml.safe_load(f2)
        return data1, data2

    def convert_diff_to_dict(self):
        """
        Purpose:    Convert the DeepDiff result to a categorized dictionary
        Modifies:   Nothing
        Returns:    Dictionary representation of the categorized DeepDiff result
        """
        differences = {}

        pattern = r"root\['(.*?)'\]"

        for change_type, changes in self.differences.items():
            if change_type in ["dictionary_item_added", "dictionary_item_removed"]:
                category = (
                    "Fields Unique to File 2"
                    if change_type == "dictionary_item_added"
                    else "Fields Unique to File 1"
                )
                if category not in differences:
                    differences[category] = []

                if isinstance(changes, dict):
                    for change in changes.keys():
                        formatted_output = (
                            re.match(pattern, change)[1]
                            if re.match(pattern, change)
                            else change
                        )
                        differences[category].append(formatted_output)
                elif isinstance(changes, list):
                    differences[category].extend(changes)
                else:
                    changes = str(changes)
                    matched_fields = re.findall(pattern, changes)
                    differences[category].extend(matched_fields)

            elif change_type == "values_changed":
                if "Values Changed" not in differences:
                    differences["Values Changed"] = {}

                for change, details in changes.items():
                    formatted_output = (
                        re.match(pattern, change)[1]
                        if re.match(pattern, change)
                        else change
                    )
                    differences["Values Changed"][
                        formatted_output
                    ] = f"{details['old_value']} -> {details['new_value']}"

            elif change_type == "type_changes":
                if "Type Changes" not in differences:
                    differences["Type Changes"] = {}

                for change, details in changes.items():
                    formatted_output = (
                        re.match(pattern, change)[1]
                        if re.match(pattern, change)
                        else change
                    )
                    differences["Type Changes"][
                        formatted_output
                    ] = f"{details['old_type']} -> {details['new_type']}"

            elif change_type == "iterable_item_added":
                if "Values Added in File 2" not in differences:
                    differences["Values Added in File 2"] = {}

                for change, details in changes.items():
                    formatted_output = (
                        re.match(pattern, change)[1]
                        if re.match(pattern, change)
                        else change
                    )
                    differences["Values Added in File 2"][formatted_output] = details

            elif change_type == "iterable_item_removed":
                if "Values Removed in File 2" not in differences:
                    differences["Values Removed in File 2"] = {}

                for change, details in changes.items():
                    formatted_output = (
                        re.match(pattern, change)[1]
                        if re.match(pattern, change)
                        else change
                    )
                    differences["Values Removed in File 2"][formatted_output] = details

        return differences
