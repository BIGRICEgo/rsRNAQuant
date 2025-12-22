
import pandas as pd
import sys
import os

def consolidate_rRNA_matrix(gse_id, base_results_path, samples, case_samples, control_samples):
    """
    Consolidates rRNA annotation results from all samples of a given GSE into a single matrix.

    Args:
        gse_id (str): The GEO accession ID (e.g., "GSE128004").
        base_results_path (str): The base path to the EVDatabase results directory.
        samples (list): List of sample names.
        case_samples (list): List of case sample names.
        control_samples (list): List of control sample names.

    Returns:
        pd.DataFrame: A DataFrame containing the consolidated matrix, or None if no samples found.
    """

    gse_stat_path = base_results_path
    output_file_path = os.path.join(gse_stat_path, f"{gse_id}_rRNA_expression_matrix.txt")

    if not samples:
        print(f"No samples found for GSE: {gse_id}")
        return None

    # Initialize a dictionary to store sequence information and RPMs
    consolidated_data = {}

    # Define the common columns for rRNA annotation
    annotation_cols = ['Sequence', 'Origin', 'Subtype', 'Fragment', 'Length', 'Start', 'End']
    
    for sample_name in samples:
        file_path = os.path.join(gse_stat_path, f"{sample_name}_rRNA_output_with_rpm.txt")
        print(f"Processing {file_path}...")
        
        if not os.path.exists(file_path):
            print(f"Warning: File not found for sample {sample_name}: {file_path}. Skipping.")
            continue

        try:
            # Read the sample's rRNA output file
            df_sample = pd.read_csv(file_path, sep='\t', header=0, comment='#') # Assuming tab-separated
            
            # Ensure required columns exist
            required_cols = annotation_cols + ['RPM']
            if not all(col in df_sample.columns for col in required_cols):
                print(f"Error: Missing required columns in {file_path}. Expected: {required_cols}. Skipping.")
                continue

            # Iterate over rows of the current sample's dataframe
            for index, row in df_sample.iterrows():
                sequence = row['Sequence']
                
                if sequence not in consolidated_data:
                    # If new sequence, add all annotation info
                    consolidated_data[sequence] = {col: row[col] for col in annotation_cols[1:]} # Exclude 'Sequence' itself
                    # Initialize RPM for all samples to 0
                    for s in samples:
                        consolidated_data[sequence][s] = 0.0
                
                # Update RPM for the current sample
                consolidated_data[sequence][sample_name] = row['RPM']
        
        except Exception as e:
            print(f"Error processing file {file_path}: {e}")
            continue

    if not consolidated_data:
        print("No data was processed. Exiting.")
        return None

    # Convert the dictionary to a DataFrame
    rows_for_df = []
    for seq, data in consolidated_data.items():
        row_dict = {'Sequence': seq}
        row_dict.update(data)
        rows_for_df.append(row_dict)

    final_df = pd.DataFrame(rows_for_df)

    # Reorder columns: annotation_cols first, then sample RPM columns
    ordered_columns = annotation_cols + samples
    # Filter out columns that might not exist if some samples were skipped
    ordered_columns = [col for col in ordered_columns if col in final_df.columns]
    final_df = final_df[ordered_columns]

    # Filter to only include samples that actually had data processed (i.e., columns exist in final_df)
    case_rpm_columns = [s for s in case_samples if s in final_df.columns]
    control_rpm_columns = [s for s in control_samples if s in final_df.columns]

    # Calculate mean and median for case samples
    if case_rpm_columns:
        final_df['mean_case'] = final_df[case_rpm_columns].mean(axis=1).round(5)
        final_df['median_case'] = final_df[case_rpm_columns].median(axis=1).round(5)
    else:
        final_df['mean_case'] = 0.0
        final_df['median_case'] = 0.0

    # Calculate mean and median for control samples
    if control_rpm_columns:
        final_df['mean_control'] = final_df[control_rpm_columns].mean(axis=1).round(5)
        final_df['median_control'] = final_df[control_rpm_columns].median(axis=1).round(5)
    else:
        final_df['mean_control'] = 0.0
        final_df['median_control'] = 0.0

    # Add new columns: count of samples where RPM is not 0, separately for case and control
    final_df['Sample_Count_case'] = final_df[case_rpm_columns].apply(lambda row: (row > 0).sum(), axis=1) if case_rpm_columns else 0
    final_df['Sample_Count_control'] = final_df[control_rpm_columns].apply(lambda row: (row > 0).sum(), axis=1) if control_rpm_columns else 0

    # Save the consolidated matrix to a text file
    try:
        final_df.to_csv(output_file_path, sep='\t', index=False)
        print(f"Successfully created consolidated matrix at: {output_file_path}")
    except Exception as e:
        print(f"Error saving output file {output_file_path}: {e}")
        return None

    return final_df

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python ExprMatrix_modified.py <GSE_ID> <BASE_RESULTS_PATH> <SAMPLES> <CASE_CONTROL>")
        sys.exit(1)

    target_gse = sys.argv[1]
    base_results_path = sys.argv[2]
    samples = sys.argv[3].split(',')
    case_control = eval(sys.argv[4]) 
    case_samples = case_control['case']
    control_samples = case_control['control']

    print(f"Consolidating rRNA matrix for GSE: {target_gse}")
    consolidated_df = consolidate_rRNA_matrix(target_gse, base_results_path, samples, case_samples, control_samples)

    if consolidated_df is not None:
        print(f"Matrix shape: {consolidated_df.shape}")
