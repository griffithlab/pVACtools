def validate_aggregated_columns(aggregated_columns, parser):
    """
    Purpose:    Makes sure the user inputs valid aggregated tsv columns
    Modifies:   Nothing
    Returns:    None
    """
    valid_aggregated_columns = [
        "Gene",
        "AA Change",
        "Num Passing Transcripts",
        "Best Peptide",
        "Best Transcript",
        "Num Passing Peptides",
        "IC50 MT",
        "IC50 WT",
        "%ile MT",
        "%ile WT",
        "RNA Expr",
        "RNA VAF",
        "DNA VAF",
        "Tier",
    ]
    for col in aggregated_columns:
        if col not in valid_aggregated_columns:
            parser.error(
                f"Invalid aggregated column '{col}' specified.\nValid columns are: {', '.join(valid_aggregated_columns)}"
            )


def validate_unaggregated_columns(unaggregated_columns, parser):
    """
    Purpose:    Makes sure the user inputs valid unaggregated tsv columns
    Modifies:   Nothing
    Returns:    None
    """
    valid_unaggregated_columns = [
        "Biotype",
        "Sub-peptide Position",
        "Median MT IC50 Score",
        "Median WT IC50 Score",
        "Median MT Percentile",
        "Median WT Percentile",
        "WT Epitope Seq",
        "Tumor DNA VAF",
        "Tumor RNA Depth",
        "Tumor RNA VAF",
        "Gene Expression",
        "BigMHC_EL WT Score",
        "BigMHC_EL MT Score",
        "BigMHC_IM WT Score",
        "BigMHC_IM MT Score",
        "MHCflurryEL Processing WT Score",
        "MHCflurryEL Processing MT Score",
        "MHCflurryEL Presentation WT Score",
        "MHCflurryEL Presentation MT Score",
        "MHCflurryEL Presentation WT Percentile",
        "MHCflurryEL Presentation MT Percentile",
        "MHCflurry WT IC50 Score",
        "MHCflurry MT IC50 Score",
        "MHCflurry WT Percentile",
        "MHCflurry MT Percentile",
        "MHCnuggetsI WT IC50 Score",
        "MHCnuggetsI MT IC50 Score",
        "MHCnuggetsI WT Percentile",
        "MHCnuggetsI MT Percentile",
        "NetMHC WT IC50 Score",
        "NetMHC MT IC50 Score",
        "NetMHC WT Percentile",
        "NetMHC MT Percentile",
        "NetMHCcons WT IC50 Score",
        "NetMHCcons MT IC50 Score",
        "NetMHCcons WT Percentile",
        "NetMHCcons MT Percentile",
        "NetMHCpan WT IC50 Score",
        "NetMHCpan MT IC50 Score",
        "NetMHCpan WT Percentile",
        "NetMHCpan MT Percentile",
        "NetMHCpanEL WT Score",
        "NetMHCpanEL MT Score",
        "NetMHCpanEL WT Percentile",
        "NetMHCpanEL MT Percentile",
        "PickPocket WT IC50 Score",
        "PickPocket MT IC50 Score",
        "PickPocket WT Percentile",
        "PickPocket MT Percentile",
        "SMM WT IC50 Score",
        "SMM MT IC50 Score",
        "SMM WT Percentile",
        "SMM MT Percentile",
        "SMMPMBEC WT IC50 Score",
        "SMMPMBEC MT IC50 Score",
        "SMMPMBEC WT Percentile",
        "SMMPMBEC MT Percentile",
        "DeepImmuno WT Score",
        "DeepImmuno MT Score",
        "Problematic Positions",
    ]
    for col in unaggregated_columns:
        if col not in valid_unaggregated_columns:
            parser.error(
                f"Invalid unaggregated column '{col}' specified.\nValid columns are: {', '.join(valid_unaggregated_columns)}"
            )


def validate_reference_match_columns(reference_match_columns, parser):
    """
    Purpose:    Makes sure the user inputs valid reference match tsv columns
    Modifies:   Nothing
    Returns:    None
    """
    valid_reference_match_columns = [
        "Peptide",
        "Hit Definition",
        "Match Window",
        "Match Sequence",
    ]
    for col in reference_match_columns:
        if col not in valid_reference_match_columns:
            parser.error(
                f"Invalid reference match column '{col}' specified.\nValid columns are: {', '.join(valid_reference_match_columns)}"
            )
