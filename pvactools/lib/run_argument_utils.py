import argparse

from pvactools.lib.run_utils import valid_tiers

def float_range(minimum, maximum):
    """Return function handle of an argument type function for
       ArgumentParser checking a float range: minimum <= arg <= maximum
         minimum - minimum acceptable argument
         maximum - maximum acceptable argument"""

    # Define the function with default arguments
    def float_range_checker(arg):
        """New Type function for argparse - a float within predefined range."""

        try:
            f = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("must be a floating point number")
        if f < minimum or f > maximum:
            raise argparse.ArgumentTypeError("must be in range [" + str(minimum) + " .. " + str(maximum)+"]")
        return f

    # Return function handle to checking function
    return float_range_checker

def aggregate_report_evaluations():
    """Return function handle of an argument type function for
       ArgumentParser checking of the aggregate report evaluation values.
       Valid values are: ['Accept', 'Reject', 'Pending', 'Review']"""

    valid_values = ['Accept', 'Reject', 'Pending', 'Review']

    def aggregate_report_evaluation_checker(arg):
        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in valid_values:
                raise argparse.ArgumentTypeError(
                    "Invalid evaluation '{}'. Valid values are: {}".format(argument, ", ".join(valid_values))
                )
        return arg_list

    return aggregate_report_evaluation_checker

def transcript_prioritization_strategy():
    """Return function handle of an argument type function for
       ArgumentParser checking of the transcript prioritization strategy
       checking that the specified criteria are in the list of: ['canonical', 'mane_select', 'tsl']"""

    # Define the function with default arguments
    def transcript_prioritization_strategy_checker(arg):
        """New Type function for argparse - a comma-separated list with predefined valid values."""

        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in ['canonical', 'mane_select', 'tsl']:
                raise argparse.ArgumentTypeError("List element must be one of 'canonical', 'mane_select', 'tsl', not {}".format(argument))
        return arg_list

    # Return function handle to checking function
    return transcript_prioritization_strategy_checker

def top_score_metric2():
    """Return function handle of an argument type function for
       ArgumentParser checking of the top score metric2
       checking that the specified criteria are in the list of: ['ic50', 'combined_percentile', 'binding_percentile', 'immunogenicity_percentile', 'presentation_percentile']"""

    # Define the function with default arguments
    def top_score_metric2_checker(arg):
        """New Type function for argparse - a comma-separated list with predefined valid values."""

        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in ['ic50', 'combined_percentile', 'binding_percentile', 'immunogenicity_percentile', 'presentation_percentile']:
                raise argparse.ArgumentTypeError("List element must be one of 'ic50', 'combined_percentile', 'binding_percentile', 'immunogenicity_percentile', 'presentation_percentile', not {}".format(argument))
        return arg_list

    # Return function handle to checking function
    return top_score_metric2_checker

def pvacsplice_anchors():
    """Return function handle of an argument type function for
       ArgumentParser checking of the pVACsplice anchors
       checking that the specified criteria are in the list of: ['A', 'D', 'NDA']"""

    # Define the function with default arguments
    def pvacsplice_anchors_checker(arg):
        """New Type function for argparse - a comma-separated list with predefined valid values."""

        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in ['A', 'D', 'NDA']:
                raise argparse.ArgumentTypeError("List element must be one of 'A', 'D', 'NDA', not {}".format(argument))
        return arg_list

    # Return function handle to checking function
    return pvacsplice_anchors_checker

def tiers(tool):
    tiers = valid_tiers(tool)
    tiers_string = ", ".join(['"{}"'.format(x) for x in tiers])
    """Return function handle of an argument type function for
       ArgumentParser checking of the pVACseq tiers
       checking that the specified criteria are in the list of: [{}]""".format(tiers_string)

    # Define the function with default arguments
    def tiers_checker(arg):
        """New Type function for argparse - a comma-separated list with predefined valid values."""

        arg_list = arg.split(",")
        for argument in arg_list:
            if argument not in tiers:
                raise argparse.ArgumentTypeError('List element must be one of {}, not {}'.format(tiers_string, argument))
        return arg_list

    # Return function handle to checking function
    return tiers_checker

def downstream_sequence_length():
    """Return function handle of an argument type function for
       ArgumentParser checking of the aggregate report evaluation values.
       Valid values are: 'full' or a positive integer"""

    def downstream_sequence_length_checker(arg):
        if arg == 'full':
            return None
        elif arg.isdigit():
            return int(arg)
        else:
            raise argparse.ArgumentTypeError("Argument needs to be a positive integer or 'full'")

    return downstream_sequence_length_checker
