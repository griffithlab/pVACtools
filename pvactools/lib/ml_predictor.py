# --------------------------------------------------------------------------------------------------
# Script: ml_predictor.py
# Description: ML prediction module for pVACtools that merges immuno pipeline output files 
#              and performs ML-based neoantigen evaluation predictions.
# --------------------------------------------------------------------------------------------------

import os
import sys
import pandas as pd
import joblib
import pickle
import numpy as np
from pathlib import Path


def merge_and_prepare_data(file1_path, file2_path, file3_path):
    """
    Merge the three input files and prepare data for ML prediction.
    
    Args:
        file1_path (str): Path to MHC Class I aggregated epitopes file
        file2_path (str): Path to MHC Class I all epitopes file  
        file3_path (str): Path to MHC Class II aggregated epitopes file
        
    Returns:
        pd.DataFrame: Prepared data ready for ML prediction
    """
    print("Reading input files...")
    
    # Read the three input files
    mhc1_agg_df = pd.read_csv(file1_path, sep='\t', na_values=["NA", "NaN", ""], keep_default_na=False)
    mhc1_allepi_df = pd.read_csv(file2_path, sep='\t', na_values=["NA", "NaN", ""], keep_default_na=False)
    mhc2_agg_df = pd.read_csv(file3_path, sep='\t', na_values=["NA", "NaN", ""], keep_default_na=False)
    
    # Rename columns in mhc1_agg_df
    mhc1_agg_df.rename(columns={
        "Best Peptide": "Best Peptide class1", 
        "Best Transcript": "Best Transcript class1", 
        "IC50 MT": "IC50 MT class1", 
        "IC50 WT": "IC50 WT class1", 
        "%ile MT": "%ile MT class1", 
        "%ile WT": "%ile WT class1"
    }, inplace=True)

    # Rename columns in mhc2_agg_df
    mhc2_agg_df.rename(columns={
        "Best Peptide": "Best Peptide class2", 
        "Best Transcript": "Best Transcript class2", 
        "IC50 MT": "IC50 MT class2", 
        "IC50 WT": "IC50 WT class2", 
        "%ile MT": "%ile MT class2", 
        "%ile WT": "%ile WT class2"
    }, inplace=True)

    # Columns to keep in mhc2_agg_df (doing this step first becuase of not, merging data would cause names such as "Index.x" and "Index.y", making column selection difficult later)
    columns_to_keep = [
        "ID",  # Always include "ID"
        *[col for col in mhc2_agg_df.columns if col.startswith("D") and "DNA VAF" not in col],  # Columns starting with "D" but without "DNA VAF"
        *[col for col in mhc2_agg_df.columns if col.endswith("class2")]  # Renamed columns ending with "class2"
    ]
    
    # Select the subset of columns
    mhc2_agg_subset = mhc2_agg_df[columns_to_keep]

    # Compare the number of rows
    if mhc1_agg_df.shape[0] != mhc2_agg_subset.shape[0]:
        print("Warning: Class 1 aggregated file DOES NOT have the same number of rows as Class 2 aggregated file.\n Extra rows in Class 1 aggregated file will have Evaluation as Pending.")

    # Merge class 1 and class 2 aggregated dataframes
    merged_df = pd.merge(mhc1_agg_df, mhc2_agg_subset, on="ID")

    # Split the "ID" column in merged_df into separate columns for matching to all epitopes file
    merged_df_split = merged_df.copy()
    merged_df_split[['Chromosome', 'Start', 'Stop', 'Reference', 'Variant']] = merged_df_split['ID'].str.split('-', expand=True)

    # Convert "Start" and "Stop" columns to integers
    merged_df_split['Start'] = merged_df_split['Start'].astype(int)
    merged_df_split['Stop'] = merged_df_split['Stop'].astype(int)
    
    # Convert Chromosome column to string in both dataframes to ensure consistent data types for merging
    merged_df_split['Chromosome'] = merged_df_split['Chromosome'].astype(str)
    mhc1_allepi_df['Chromosome'] = mhc1_allepi_df['Chromosome'].astype(str)

    # Perform an inner join with mhc1_allepi_df on the specified columns
    merged_all = pd.merge(
        merged_df_split,
        mhc1_allepi_df,
        how='inner',
        left_on=[
            'Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 
            'Best Transcript class1', 'Best Peptide class1', 'Allele'
        ],
        right_on=[
            'Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 
            'Transcript', 'MT Epitope Seq', 'HLA Allele'
        ]
    )

    # Drop the redundant columns created from splitting "ID"
    merged_all.drop(columns=["Chromosome", "Start", "Stop", "Reference", "Variant", "Transcript", "MT Epitope Seq", "HLA Allele"], inplace=True)

    # Transformations
    # NOTE: Pos may take form "#-#", so we need to extract the first integer
    if merged_all["Pos"].dtype == "object":
        merged_all["Pos"] = merged_all["Pos"].str.extract(r"^(\d+)").astype("Int64")  # Extract the first integer from Pos
    # Convert Pos to float64
    merged_all['Pos'] = merged_all['Pos'].astype(float)

    # NOTE: Prob Pos may take form "#,#", so we need to extract the first integer
    merged_all["Prob Pos"] = (
        merged_all["Prob Pos"]
        .fillna("0")  # Replace NaN with "0"
        .astype(str)  # Ensure all values are strings
        .replace("None", "0")  # Replace "None" with "0"
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

    # Select columns to keep from merged_all
    columns_keep = [
        "ID",
        "Evaluation",
        "TSL",
        "Pos",
        "Prob Pos",
        "Num Passing Peptides",
        "RNA Expr",
        "RNA VAF",
        "Allele Expr",
        "RNA Depth",
        "DNA VAF",
        "Ref Match",
        "Biotype",
        "Variant Type",
        "Peptide Length",
        "Best MT IC50 Score",
        "Corresponding WT IC50 Score",
        "Corresponding Fold Change",
        "Best MT Percentile",
        "Corresponding WT Percentile",
        "Median Fold Change",
        "cysteine_count"
    ]

    # Select the final columns
    final_columns = columns_keep + [
        col for col in merged_all.columns if col.startswith(("IC50", "%ile", "MHCflurry", "MHCnuggetsI", "NetMHC", "PickPocket", "SMM"))
    ] + ['Prob match', 'Gene of Interest']

    merged_all = merged_all[final_columns]

    # Unify the column names
    merged_all.rename(columns={
        "NetMHCpanEL WT Score": "NetMHCpanEL WT IC50 Score",
        "NetMHCpanEL MT Score": "NetMHCpanEL MT IC50 Score",
        "MHCflurry WT Score": "MHCflurry WT IC50 Score",
        "MHCflurry MT Score": "MHCflurry MT IC50 Score",
        "MHCnuggetsI WT Score": "MHCnuggetsI WT IC50 Score",
        "MHCnuggetsI MT Score": "MHCnuggetsI MT IC50 Score",
        "NetMHC WT Score": "NetMHC WT IC50 Score",
        "NetMHC MT Score": "NetMHC MT IC50 Score",
        "NetMHCcons WT Score": "NetMHCcons WT IC50 Score",
        "NetMHCcons MT Score": "NetMHCcons MT IC50 Score",
        "NetMHCpan WT Score": "NetMHCpan WT IC50 Score",
        "NetMHCpan MT Score": "NetMHCpan MT IC50 Score",
        "PickPocket WT Score": "PickPocket WT IC50 Score",
        "PickPocket MT Score": "PickPocket MT IC50 Score",
        "SMM WT Score": "SMM WT IC50 Score",
        "SMM MT Score": "SMM MT IC50 Score",
        "SMMPMBEC WT Score": "SMMPMBEC WT IC50 Score",
        "SMMPMBEC MT Score": "SMMPMBEC MT IC50 Score"
    }, inplace=True)

    return merged_all

def _resolve_artifact_paths(model_artifacts_path):
    """
    Resolve standard artifact file paths inside the provided artifacts directory.

    Expected filenames (numpy126 versions):
      - rf_downsample_model_numpy126.pkl
      - trained_imputer_numpy126.joblib
      - label_encoders_numpy126.pkl
    """
    artifacts_dir = Path(model_artifacts_path)
    model_path = artifacts_dir / 'rf_downsample_model_numpy126.pkl'
    imputer_path = artifacts_dir / 'trained_imputer_numpy126.joblib'
    encoders_path = artifacts_dir / 'label_encoders_numpy126.pkl'
    return str(model_path), str(imputer_path), str(encoders_path)

def clean_and_impute_data(merged_all, model_artifacts_path):
    """
    Clean and impute the merged data for ML prediction.
    
    Args:
        merged_all (pd.DataFrame): Merged data from merge_and_prepare_data
        imputer_path (str): Path to the trained imputer
        label_encoders_path (str): Path to the label encoders
        
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

def make_ml_predictions(post_imputed_data, model_artifacts_path, threshold=0.55):
    """
    Make ML predictions using the trained model.
    
    Args:
        post_imputed_data (pd.DataFrame): Cleaned and imputed data
        model_path (str): Path to the trained ML model
        threshold (float): Threshold for Accept predictions (default: 0.55)
        
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
            post_imputed_data['Accept_pred_prob'] >= threshold, "Accept",
            np.where(
                post_imputed_data['Accept_pred_prob'] > 0.30, "Review", "Reject" # keep "Review" as "Review", change to Pending in a later step
            )
        )
    )
    
    return post_imputed_data

def create_final_output(post_imputed_data, original_agg_file_path, output_dir, sample_name):
    """
    Create the final output files with ML predictions.
    
    Args:
        post_imputed_data (pd.DataFrame): Data with ML predictions
        original_agg_file_path (str): Path to original aggregated file
        output_dir (str): Output directory for results
        sample_name (str): Sample name for output file
        
    Returns:
        tuple: Paths to the two final output files (comments format, new format)
    """
    print("Creating final output files...")
    
    # Read original aggregated file
    mhc1_agg_df_itb = pd.read_csv(original_agg_file_path, sep="\t", dtype=str)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Version 1: put the prediction and score in the Comments column
    final_df = mhc1_agg_df_itb.drop(columns=['Evaluation']).merge(
        post_imputed_data[['ID', 'Evaluation_pred', 'Accept_pred_prob']], on="ID", how="left"
    ).rename(columns={'Evaluation_pred': 'Evaluation'})
    final_df['Comments'] = "Probability of Accept: " + final_df['Accept_pred_prob'].round(3).astype(str)
    final_df = final_df.drop(columns=['Accept_pred_prob'])
    final_df.loc[final_df['Evaluation'] == "Pending", 'Comments'] = "Unable to make prediction with ML model"
    # Modify 'Review' to 'Pending' in Evaluation column
    final_df.loc[final_df['Evaluation'] == 'Review', 'Evaluation'] = 'Pending'
    
    # Version 2: put the prediction and score in an extra column, instead of the Comments column
    final_df2 = mhc1_agg_df_itb.drop(columns=['Evaluation']).merge(
        post_imputed_data[['ID', 'Evaluation_pred', 'Accept_pred_prob']], on="ID", how="left"
    ).rename(columns={'Evaluation_pred': 'Evaluation'})

    # Conditionally set ML Prediction (score) based on Evaluation_pred
    # NOTE: sometimes the ML model will not be able to make a prediction due to missing data, in this case, the Evaluation_pred will be "Pending"
    # This may be due to the class 1 and class 2 files not having the same number of rows.
    final_df2['ML Prediction (score)'] = final_df2.apply(
        lambda row: "Unable to make prediction with ML model due to missing data" if row['Evaluation'] == "Pending" 
        else row['Evaluation'] + " (" + str(round(row['Accept_pred_prob'], 2)) + ")", 
        axis=1
    )
    # Modify 'Review' to 'Pending' in Evaluation column
    final_df2.loc[final_df2['Evaluation'] == 'Review', 'Evaluation'] = 'Pending'
    final_df2 = final_df2.drop(columns=['Accept_pred_prob'])
    
    # Save both output files
    output_file1 = os.path.join(output_dir, f"{sample_name}_predict_pvacview.tsv")
    output_file2 = os.path.join(output_dir, f"{sample_name}_predict_pvacview_new_format.tsv")
    
    final_df.to_csv(output_file1, sep="\t", index=False, na_rep="NA")
    final_df2.to_csv(output_file2, sep="\t", index=False, na_rep="NA")
    
    print(f"ML predictions saved to: {output_file1} and {output_file2}")
    return output_file1, output_file2

def run_ml_predictions(file1_path, file2_path, file3_path, model_artifacts_path, output_dir, sample_name, threshold=0.55):
    """
    Main function to run ML predictions on pVACtools output files.
    
    Args:
        file1_path (str): Path to MHC Class I aggregated epitopes file
        file2_path (str): Path to MHC Class I all epitopes file
        file3_path (str): Path to MHC Class II aggregated epitopes file
        model_artifacts_path (str): Path to directory containing model artifacts
        output_dir (str): Output directory for results
        sample_name (str): Sample name for output file
        threshold (float): Threshold for Accept predictions (default: 0.55)
        
    Returns:
        tuple: Paths to the two final output files (comments format, new format)
    """
    print("Starting ML prediction pipeline...")
    
    try:
        # Step 1: Merge and prepare data
        merged_all = merge_and_prepare_data(file1_path, file2_path, file3_path)
        
        # Step 2: Clean and impute data
        post_imputed_data = clean_and_impute_data(merged_all, model_artifacts_path)
        
        # Step 3: Make ML predictions
        post_imputed_data = make_ml_predictions(post_imputed_data, model_artifacts_path, threshold)
        
        # Step 4: Create final output
        output_files = create_final_output(post_imputed_data, file1_path, output_dir, sample_name)
        
        print("ML prediction pipeline completed successfully!")
        return output_files
        
    except Exception as e:
        print(f"Error in ML prediction pipeline: {str(e)}")
        raise

if __name__ == "__main__":
    # This section is for testing the script independently
    import argparse
    
    parser = argparse.ArgumentParser(description="ML prediction for pVACtools output")
    parser.add_argument("--file1", required=True, help="MHC Class I aggregated epitopes file")
    parser.add_argument("--file2", required=True, help="MHC Class I all epitopes file")
    parser.add_argument("--file3", required=True, help="MHC Class II all epitopes file")
    parser.add_argument("--artifacts", required=True, help="Path to directory containing model artifacts")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--threshold", type=float, default=0.55, help="Prediction threshold")
    
    args = parser.parse_args()
    
    run_ml_predictions(
        args.file1, args.file2, args.file3,
        args.artifacts,
        args.output, args.sample, args.threshold
    )
