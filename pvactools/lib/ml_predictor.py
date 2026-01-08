# --------------------------------------------------------------------------------------------------
# Script: ml_predictor.py
# Description: ML prediction module for pVACtools that merges immuno pipeline output files 
#              and performs ML-based neoantigen evaluation predictions.
# --------------------------------------------------------------------------------------------------

import argparse
import os
import pandas as pd
import joblib
import pickle
import numpy as np
from pathlib import Path


def merge_and_prepare_data(class1_aggregated_path, class1_all_epitopes_path, class2_aggregated_path):
    """
    Merge the three input files and prepare data for ML prediction.
    
    Args:
        class1_aggregated_path (str): Path to MHC Class I aggregated epitopes file
        class1_all_epitopes_path (str): Path to MHC Class I all epitopes file  
        class2_aggregated_path (str): Path to MHC Class II aggregated epitopes file
        
    Returns:
        pd.DataFrame: Prepared data ready for ML prediction
    """
    print("Reading input files...")
    
    # --- Columns to read from input files ---
    mhc1_agg_columns = [
        'ID', 'Index', 'Best Peptide', 'IC50 MT', 'IC50 WT', '%ile MT', '%ile WT',
        'Allele', 'Allele Expr', 'DNA VAF', 'TSL', 'Num Passing Peptides', 
        'Pos', 'Prob Pos', 'RNA Depth', 'RNA Expr', 'RNA VAF', 
        'Ref Match', 'Evaluation'
    ]
    mhc2_agg_columns = [
        'ID', '%ile MT', '%ile WT', 'IC50 MT', 'IC50 WT'
    ]
    mhc1_allepi_columns = [
        'Index', 'MT Epitope Seq', 'HLA Allele', 'Gene of Interest', 
        'Corresponding Fold Change', 'Variant Type', 'Best MT IC50 Score', 
        'Best MT Percentile', 'Corresponding WT IC50 Score', 
        'Corresponding WT Percentile', 'Biotype', 'cysteine_count'
    ]
    predictor_columns = [
        'MHCflurry MT IC50 Score',
        'MHCflurry MT Percentile',
        'MHCflurry WT IC50 Score',
        'MHCflurry WT Percentile',
        'MHCflurryEL Presentation MT Percentile',
        'MHCflurryEL Presentation MT Score',
        'MHCflurryEL Presentation WT Percentile',
        'MHCflurryEL Presentation WT Score',
        'MHCflurryEL Processing MT Score',
        'MHCflurryEL Processing WT Score',
        'MHCnuggetsI MT IC50 Score',
        'MHCnuggetsI MT Percentile',
        'MHCnuggetsI WT IC50 Score',
        'MHCnuggetsI WT Percentile',
        'Median Fold Change',
        'NetMHC MT IC50 Score',
        'NetMHC MT Percentile',
        'NetMHC WT IC50 Score',
        'NetMHC WT Percentile',
        'NetMHCcons MT IC50 Score',
        'NetMHCcons MT Percentile',
        'NetMHCcons WT IC50 Score',
        'NetMHCcons WT Percentile',
        'NetMHCpan MT IC50 Score',
        'NetMHCpan MT Percentile',
        'NetMHCpan WT IC50 Score',
        'NetMHCpan WT Percentile',
        'NetMHCpanEL MT Presentation Score',
        'NetMHCpanEL MT Percentile',
        'NetMHCpanEL WT Presentation Score',
        'NetMHCpanEL WT Percentile',
        'Peptide Length',
        'PickPocket MT IC50 Score',
        'PickPocket MT Percentile',
        'PickPocket WT IC50 Score',
        'PickPocket WT Percentile',
        'SMM MT IC50 Score',
        'SMM MT Percentile',
        'SMM WT IC50 Score',
        'SMM WT Percentile',
        'SMMPMBEC MT IC50 Score',
        'SMMPMBEC MT Percentile',
        'SMMPMBEC WT IC50 Score',
        'SMMPMBEC WT Percentile'
    ]
    # --- Read the three input files ---
    mhc1_agg_df = pd.read_csv(class1_aggregated_path, sep='\t', na_values=["NA", "NaN", ""], keep_default_na=False, usecols=mhc1_agg_columns)
    mhc1_allepi_df = pd.read_csv(class1_all_epitopes_path, sep='\t', na_values=["NA", "NaN", ""], keep_default_na=False, usecols=mhc1_allepi_columns+predictor_columns)
    mhc2_agg_df = pd.read_csv(class2_aggregated_path, sep='\t', na_values=["NA", "NaN", ""], keep_default_na=False, usecols=mhc2_agg_columns)
    
    # Rename columns in mhc1_agg_df
    mhc1_agg_df.rename(columns={ 
        "IC50 MT": "IC50 MT class1", 
        "IC50 WT": "IC50 WT class1", 
        "%ile MT": "%ile MT class1", 
        "%ile WT": "%ile WT class1"
    }, inplace=True)

    # Rename columns in mhc2_agg_df
    mhc2_agg_df.rename(columns={
        "IC50 MT": "IC50 MT class2", 
        "IC50 WT": "IC50 WT class2", 
        "%ile MT": "%ile MT class2", 
        "%ile WT": "%ile WT class2"
    }, inplace=True)

    # Compare the number of rows
    if mhc1_agg_df.shape[0] != mhc2_agg_df.shape[0]:
        print("Warning: Class 1 aggregated file DOES NOT have the same number of rows as Class 2 aggregated file.\n May cause \"NA\" in ML predictions.")

    # --- Merge class 1 and class 2 aggregated dataframes ---
    merged_df = pd.merge(mhc1_agg_df, mhc2_agg_df, on="ID")

    # Perform an inner join with mhc1_allepi_df on the specified columns ---
    merged_all = pd.merge(
        merged_df,
        mhc1_allepi_df,
        how='inner',
        left_on=['Index', 'Best Peptide', 'Allele'],
        right_on=['Index', 'MT Epitope Seq', 'HLA Allele']
    )

    # --- Transformations ---
    # NOTE: Pos may take form "#-#" or "#,#", so we need to extract the first integer
    # NOTE: the "#-#" format is from older version of pVACtools that I encountered during training the model, it might not be needed anymore
    if merged_all["Pos"].dtype == "object":
        # Extract the first integer (handles formats like "5", "5-10", "5,10", etc.)
        merged_all["Pos"] = merged_all["Pos"].astype(str).str.extract(r"^(\d+)").astype("Int64")
        # Convert Pos to float64
        merged_all['Pos'] = merged_all['Pos'].astype(float)
    else:
        # Convert Pos to float64
        merged_all['Pos'] = merged_all['Pos'].astype(float)

    # NOTE: Prob Pos may take form "#,#", so we need to extract the first integer
    # NOTE: these formats might also be from older version of pVACtools 
    merged_all["Prob Pos"] = (
        merged_all["Prob Pos"]
        .fillna("0")  # Replace NaN with "0"
        .astype(str)  # Ensure all values are strings
        .str.split(",")  # Split by commas
        .apply(lambda x: int(float(x[0])) if x[0].replace('.', '', 1).isdigit() else 0)  # Handle floats and integers
    )

    # Create a new column `Prob.match` based on whether the integer(s) in `Prob.Pos` match the value in the `Pos` column
    merged_all["Prob match"] = merged_all.apply(
        lambda row: "True" if pd.notna(row["Pos"]) and row["Pos"] in [int(x) for x in str(row["Prob Pos"]).split(",") if x.replace('.', '', 1).isdigit()] else "False",
        axis=1
    )
    
    # If they are currently strings "True"/"False", convert them:
    merged_all['Prob match'] = merged_all['Prob match'].map({'True': True, 'False': False}).astype(bool)
    merged_all['Gene of Interest'] = merged_all['Gene of Interest'].map({'True': True, 'False': False}).astype(bool)
    
    # Replace NA in TSL with 6, change into int (updated to match merge_data_predict_ml.py)
    merged_all["TSL"] = merged_all["TSL"].fillna(6).astype(int)

    # --- Rename columns for pvactools v7, otherwise I will need to retrain the model to use new column names (will do this in the future) ---
    rename_map_v7 = {
        "NetMHCpanEL MT Presentation Score": "NetMHCpanEL MT IC50 Score",
        "NetMHCpanEL WT Presentation Score": "NetMHCpanEL WT IC50 Score"
    }
    merged_all = merged_all.rename(columns=rename_map_v7)

    # --- Drop the redundant columns ---
    merged_all = merged_all.drop(columns=["Index", "MT Epitope Seq", "HLA Allele", 'Best Peptide', 'Allele'])

    return merged_all

def _get_default_artifacts_dir():
    """
    Get the path to the default ML model artifacts directory included in the pvactools package.
    
    Returns:
        str: Path to the artifacts directory
    """
    base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
    artifacts_dir = os.path.join(base_dir, 'supporting_files', 'ml_model_artifacts')
    return artifacts_dir

def _resolve_artifact_paths(model_artifacts_path=None):
    """
    Resolve standard artifact file paths inside the provided artifacts directory.
    If no path is provided, uses the default artifacts directory from the package.

    Expected filenames (numpy126 versions):
      - rf_downsample_model_numpy126.pkl
      - trained_imputer_numpy126.joblib
      - label_encoders_numpy126.pkl
    
    Args:
        model_artifacts_path (str, optional): Path to artifacts directory. If None, uses default package location.
    """
    if model_artifacts_path is None:
        artifacts_dir = _get_default_artifacts_dir()
    else:
        artifacts_dir = model_artifacts_path
    artifacts_dir = Path(artifacts_dir)
    model_path = artifacts_dir / 'rf_downsample_model_numpy126.pkl'
    imputer_path = artifacts_dir / 'trained_imputer_numpy126.joblib'
    encoders_path = artifacts_dir / 'label_encoders_numpy126.pkl'
    return str(model_path), str(imputer_path), str(encoders_path)

def clean_and_impute_data(merged_all, model_artifacts_path):
    """
    Clean and impute the merged data for ML prediction.
    
    Args:
        merged_all (pd.DataFrame): Merged data from merge_and_prepare_data
        model_artifacts_path (str): Path to the model artifacts directory
        
    Returns:
        pd.DataFrame: Cleaned and imputed data ready for prediction
    """
    print("Loading ML model artifacts...")
    
    # Load trained imputer and label encoders
    _model_path, imputer_path, label_encoders_path = _resolve_artifact_paths(model_artifacts_path)
    imputer = joblib.load(imputer_path)
    with open(label_encoders_path, "rb") as f:
        label_encoders = pickle.load(f)

    # Prepare columns for imputation
    exclude_columns = ["ID", "Evaluation"]  # Add any other columns you want to exclude
    columns_to_impute = merged_all.columns.difference(exclude_columns)

    excluded_data = merged_all[exclude_columns].copy()
    data_to_impute = merged_all[columns_to_impute].copy()

    # Apply label encoding to categorical columns
    categorical_columns = data_to_impute.select_dtypes(include=['category', 'object']).columns
    for col in categorical_columns:
        if col in label_encoders:
            le = label_encoders[col]
            # Handle unseen categories by using the most common category
            data_to_impute.loc[:, col] = data_to_impute[col].map(
                lambda x: le.transform([x])[0] if x in le.classes_ else le.transform([le.classes_[0]])[0]
            )

    # Impute missing values
    imputed_data = imputer.transform(data_to_impute)
    imputed_data = pd.DataFrame(imputed_data, columns=columns_to_impute)

    # Combine imputed and excluded columns
    post_imputed_data = pd.concat([excluded_data.reset_index(drop=True), imputed_data.reset_index(drop=True)], axis=1)
    
    return post_imputed_data

def make_ml_predictions(post_imputed_data, model_artifacts_path, threshold_accept=0.55, threshold_reject=0.30):
    """
    Make ML predictions using the trained model.
    
    Args:
        post_imputed_data (pd.DataFrame): Cleaned and imputed data
        model_path (str): Path to the trained ML model
        threshold_accept (float): Threshold for Accept predictions (default: 0.55)
        threshold_reject (float): Threshold for Reject predictions (default: 0.30)
        
    Returns:
        pd.DataFrame: Data with ML predictions added
    """
    print("Loading ML model and making predictions...")
    
    # Load the trained Random Forest model
    model_path, _imputer_path, _encoders_path = _resolve_artifact_paths(model_artifacts_path)
    rf_model = joblib.load(model_path)

    # Prepare data for prediction
    # If 'ID' or 'Evaluation' are in post_imputed_data, drop them for prediction
    predict_cols = [col for col in post_imputed_data.columns if col not in ['ID', 'Evaluation', 'patient_id']]

    # Make prediction with rf_model
    rf_pred = rf_model.predict_proba(post_imputed_data[predict_cols])[:, 1]

    # Add predictions to DataFrame
    post_imputed_data['Accept_pred_prob'] = rf_pred
    post_imputed_data['Evaluation_pred'] = np.where(
        post_imputed_data['Accept_pred_prob'].isna(), "Pending", # Set to "Pending" if the ML model is unable to make a prediction due to missing data
        np.where(
            post_imputed_data['Accept_pred_prob'] >= threshold_accept, "Accept",
            np.where(
                post_imputed_data['Accept_pred_prob'] > threshold_reject, "Review", "Reject" # keep "Review" as "Review", change to Pending in a later step
            )
        )
    )
    
    return post_imputed_data

def create_final_output(post_imputed_data, original_agg_file_path, output_dir, sample_name):
    """
    Create the final output file with ML predictions.
    
    Args:
        post_imputed_data (pd.DataFrame): Data with ML predictions
        original_agg_file_path (str): Path to original aggregated file
        output_dir (str): Output directory for results
        sample_name (str): Sample name for output file
        
    Returns:
        str: Path to the final output file
    """
    print("Creating final output file...")
    
    # 1) Read original aggregated file as pure text, preserving None/empty values
    orig_df = pd.read_csv(original_agg_file_path, sep="\t", dtype=str, keep_default_na=False, na_values=[])
    original_columns = list(orig_df.columns)
    
    # We know we will replace 'Evaluation', so treat everything else as "must preserve"
    cols_to_preserve = [c for c in original_columns if c != "Evaluation"]
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # 2) Merge in ML predictions
    #    Work from a copy that drops the old Evaluation column
    base_df = orig_df.drop(columns=["Evaluation"])
    merged_df = base_df.merge(
        post_imputed_data[['ID', 'Evaluation_pred', 'Accept_pred_prob']], 
        on="ID", 
        how="left"
    ).rename(columns={'Evaluation_pred': 'Evaluation'})
    # Fill missing evaluations from unmatched IDs as Pending
    merged_df['Evaluation'] = merged_df['Evaluation'].fillna('Pending')

    # Conditionally set ML Prediction (score) based on Evaluation_pred
    # NOTE: sometimes the ML model will not be able to make a prediction due to missing data, in this case, the Evaluation_pred will be "Pending"
    # This may be due to the class 1 and class 2 files not having the same number of rows.
    merged_df['ML Prediction (score)'] = merged_df.apply(
        lambda row: "NA" if (pd.isna(row['Evaluation']) or row['Evaluation'] == "Pending" or pd.isna(row['Accept_pred_prob']))
        else str(row['Evaluation']) + " (" + str(round(row['Accept_pred_prob'], 2)) + ")", 
        axis=1
    )
    # Modify 'Review' to 'Pending' in Evaluation column
    merged_df.loc[merged_df['Evaluation'] == 'Review', 'Evaluation'] = 'Pending'
    final_df = merged_df.drop(columns=['Accept_pred_prob'])
    
    final_df[cols_to_preserve] = orig_df[cols_to_preserve]
    
    # Save output file
    output_file = os.path.join(output_dir, f"{sample_name}_predict_pvacview.tsv")

    final_df.to_csv(output_file, float_format='%.3f', sep="\t", index=False, na_rep="NA")
    
    print(f"ML predictions saved to: {output_file}")
    return output_file

def run_ml_predictions(class1_aggregated_path, class1_all_epitopes_path, class2_aggregated_path, model_artifacts_path=None, output_dir=None, sample_name=None, threshold_accept=0.55, threshold_reject=0.30):
    """
    Main function to run ML predictions on pVACtools output files.
    
    Args:
        class1_aggregated_path (str): Path to MHC Class I aggregated epitopes file
        class1_all_epitopes_path (str): Path to MHC Class I all epitopes file
        class2_aggregated_path (str): Path to MHC Class II aggregated epitopes file
        model_artifacts_path (str, optional): Path to directory containing model artifacts. 
            If None, uses the default artifacts directory from the pvactools package.
        output_dir (str, optional): Output directory for results
        sample_name (str, optional): Sample name for output file
        threshold_accept (float): Threshold for Accept predictions (default: 0.55)
        threshold_reject (float): Threshold for Reject predictions (default: 0.30)
        
    Returns:
        str: Path to the final output file
    """
    # Validate threshold parameters
    if threshold_reject > threshold_accept:
        raise ValueError(
            f"threshold_reject ({threshold_reject}) must be less than or equal to "
            f"threshold_accept ({threshold_accept}). Please adjust your thresholds."
        )
    
    print("Starting ML prediction pipeline...")
    
    try:
        # Step 1: Merge and prepare data
        merged_all = merge_and_prepare_data(class1_aggregated_path, class1_all_epitopes_path, class2_aggregated_path)
        
        # Step 2: Clean and impute data
        post_imputed_data = clean_and_impute_data(merged_all, model_artifacts_path)
        
        # Step 3: Make ML predictions
        post_imputed_data = make_ml_predictions(post_imputed_data, model_artifacts_path, threshold_accept, threshold_reject)
        
        # Step 4: Create final output
        output_files = create_final_output(post_imputed_data, class1_aggregated_path, output_dir, sample_name)
        
        print("ML prediction pipeline completed successfully!")
        return output_files
        
    except Exception as e:
        print(f"Error in ML prediction pipeline: {str(e)}")
        raise


def define_add_ml_predictions_parser(tool='pvacseq'):
    """
    Create and return an argparse parser for the add_ml_predictions command.
    """
    parser = argparse.ArgumentParser(
        f"{tool} add_ml_predictions",
        description="Add ML-based neoantigen evaluation predictions to existing pVACseq output files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "class1_aggregated",
        help="Path to the MHC Class I aggregated epitopes TSV."
    )
    parser.add_argument(
        "class1_all_epitopes",
        help="Path to the MHC Class I all epitopes TSV."
    )
    parser.add_argument(
        "class2_aggregated",
        help="Path to the MHC Class II aggregated epitopes TSV."
    )
    parser.add_argument(
        "--artifacts_path",
        dest="artifacts_path",
        help="Optional path to a directory containing ML model artifacts. Defaults to the package-provided artifacts."
    )
    parser.add_argument(
        "output_dir",
        help="Directory where the ML prediction TSV files should be written."
    )
    parser.add_argument(
        "sample_name",
        help="Sample name prefix to use for the output files."
    )
    parser.add_argument(
        "--threshold_accept",
        type=float,
        default=0.55,
        help="Prediction threshold for Accept predictions (default: 0.55)."
    )
    parser.add_argument(
        "--threshold_reject",
        type=float,
        default=0.30,
        help="Prediction threshold for Reject predictions (default: 0.30)."
    )
    return parser


if __name__ == "__main__":
    parser = define_add_ml_predictions_parser(tool='ml_predictor')
    args = parser.parse_args()

    run_ml_predictions(
        args.class1_aggregated,
        args.class1_all_epitopes,
        args.class2_aggregated,
        args.artifacts_path,
        args.output_dir,
        args.sample_name,
        args.threshold_accept,
        args.threshold_reject
    )
