import pandas as pd
import sys
import os

def calculate_rpm(GSE, sample_name, output_dir):
    counts_file = os.path.join(output_dir, f"{sample_name}_rRNA_output.txt")
    
    # Read counts data
    df_counts = pd.read_csv(counts_file, sep='\t')

    summary_file = os.path.join(output_dir, f"{sample_name}_summary.txt")
    # Read summary file to get clean_reads
    with open(summary_file, 'r') as f:
        for line in f:
            if line.startswith('Clean_Reads'):
                clean_reads = float(line.strip().split('\t')[2])
                break
    
    # Calculate RPM
    df_counts['RPM'] = ((df_counts['Counts'] / clean_reads) * 1_000_000).round(5)
    df_counts['Sample'] = f"{sample_name}"
    df_counts['GSE'] = f"{GSE}"

    # Save df_counts with RPM to a new file
    output_counts_file = os.path.join(output_dir, f"{sample_name}_rRNA_output_with_rpm.txt")
    print(f"Successfully generated {output_counts_file} .")
    df_counts.to_csv(output_counts_file, sep='\t', index=False, header=True)

    # Output Length Distribution
    subtype_length_counts = df_counts.groupby(['Subtype', 'Length'])['Counts'].sum().reset_index()
    subtype_length_counts.rename(columns={'Counts': 'Total_Counts'}, inplace=True)
    subtype_length_counts['RPM'] = (subtype_length_counts['Total_Counts'] / clean_reads) * 1_000_000

    output_file = os.path.join(output_dir, f"{sample_name}_rRNA_length_distribution.txt")
    print(f"Successfully generated {output_file} .")

    # Save to output file: Only include 'Subtype', 'Length', and 'RPM' columns
    subtype_length_counts[['Subtype', 'Length', 'RPM']].to_csv(output_file, sep='\t', index=False, header=True)
   
    print(f"Processing for sample {sample_name} completed. Results in {output_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python RPM.py <GSE> <sample_name> <output_dir>")
        sys.exit(1)
    
    GSE = sys.argv[1]
    sample_name = sys.argv[2]
    output_dir = sys.argv[3]
    calculate_rpm(GSE, sample_name, output_dir)
