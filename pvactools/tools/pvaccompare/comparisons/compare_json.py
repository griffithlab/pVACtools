import json


class CompareJSON:
    def __init__(self, input_file1, input_file2):
        self.input_file1 = input_file1
        self.input_file2 = input_file2
        self.json1, self.json2 = self.load_files()
        self.differences = {}

    def load_files(self):
        """
        Purpose:    Load the two json metrics files into dictionaries
        Modifies:   Nothing
        Returns:    Two dictionaries corresponding to the two input files
        """
        with open(self.input_file1) as f1, open(self.input_file2) as f2:
            json1 = json.load(f1)
            json2 = json.load(f2)
        return json1, json2

    @staticmethod
    def filter_chr_keys(data):
        """
        Purpose:    Filter the data to remove everything except the inputs
        Modifies:   Nothing
        Returns:    Filtered data
        """
        if isinstance(data, dict):
            return {
                k: CompareJSON.filter_chr_keys(v)
                for k, v in data.items()
                if not k.startswith("chr")
            }
        elif isinstance(data, list):
            return [CompareJSON.filter_chr_keys(item) for item in data]
        else:
            return data

    def compare_metric_data(self):
        """
        Purpose:    Get and store all of the differences
        Modifies:   self.differences
        Returns:    None
        """
        filtered_data1 = self.filter_chr_keys(self.json1)
        filtered_data2 = self.filter_chr_keys(self.json2)

        differences = {
            "Shared Fields": sorted(
                k for k in filtered_data1.keys() if k in filtered_data2
            ),
            "Fields Unique to File 1": {
                k: v
                for k, v in sorted(filtered_data1.items())
                if k not in filtered_data2
            },
            "Fields Unique to File 2": {
                k: v
                for k, v in sorted(filtered_data2.items())
                if k not in filtered_data1
            },
            "Values Changed": {
                k: f"{sorted(filtered_data1[k])} -> {sorted(filtered_data2[k])}"
                if isinstance(filtered_data1[k], list) and isinstance(filtered_data2[k], list)
                else f"{filtered_data1[k]} -> {filtered_data2[k]}"
                for k in sorted(filtered_data1)
                if k in filtered_data2 and (
                    sorted(filtered_data1[k]) != sorted(filtered_data2[k])
                    if isinstance(filtered_data1[k], list) and isinstance(filtered_data2[k], list)
                    else filtered_data1[k] != filtered_data2[k]
                )
            },
        }
        if any(key != "Shared Fields" and differences[key] for key in differences):
            self.differences = differences
