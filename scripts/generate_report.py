import sys
import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
import yaml 

# Get input parameters
GSE = sys.argv[1]
configfile = sys.argv[2]
results_dir = sys.argv[3]
scripts_dir = sys.argv[4]
output_html_path = sys.argv[5]

# Define template environment
template_env = Environment(loader=FileSystemLoader(scripts_dir))
template = template_env.get_template('template_report.html')

# Get configuration file path
config_file = configfile

# Get sample count, case count, control count
with open(config_file, 'r') as f:
    config_data = yaml.safe_load(f)
sample_count = len(config_data['samples'])
case_count = len(config_data['groups'].get('case', []))
control_count = len(config_data['groups'].get('control', []))

# Get current report time
report_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# Define file paths
#alignment_summary_path = os.path.join(results_dir, GSE, "stats", f"{GSE}_alignment_summary.csv")
expression_matrix_path = os.path.join(results_dir, GSE, "stats", f"{GSE}_rRNA_expression_matrix.txt")
count_matrix_path = os.path.join(results_dir, GSE, "stats", f"{GSE}_rRNA_count_matrix.txt")
da_summary_path = os.path.join(results_dir, GSE, "plots", f"{GSE}_DA_summary.txt")
biomarker_table_path = os.path.join(results_dir, GSE, "plots", f"{GSE}_lasso_selected_features_rpm.txt")

# Read Differential Analysis Summary File
try:
    da_summary = pd.read_csv(da_summary_path, delimiter='\t', header=None)
    up_value = da_summary.iloc[3, 1]
    down_value = da_summary.iloc[4, 1]
except FileNotFoundError:
    up_value = "Not Enough Samples"
    down_value = "Not Enough Samples"

# Read Biomarker Table File
try:
    biomarker_table = open(biomarker_table_path, 'r').read()
except FileNotFoundError:
    biomarker_table = "Not Enough Samples"

# Render HTML template
output_html = template.render(
    report_time=report_time,
    config_file=config_file,
    sample_count=sample_count,
    case_count=case_count,
    control_count=control_count,
    results_dir=results_dir,
    GSE=GSE,
    up_value=up_value,
    down_value=down_value,
    biomarker_table=biomarker_table
)

# Write into HTML file
with open(output_html_path, 'w') as f:
    f.write(output_html)
