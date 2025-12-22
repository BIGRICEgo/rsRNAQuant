#!/bin/bash
# Usage: bash re_logs.sh <GSE> <results_dir> <samples>

set -euo pipefail

GSE="$1"
RESULTS_DIR="$2"
SAMPLES="$3"
NCRNA_FLAG="$4"
FILES_FLAG="$5"

SAMPLES="${SAMPLES//,/ }"
LOG_DIR="${RESULTS_DIR}/${GSE}/logs"

echo "[INFO] Start merging RNA logs for GSE: ${GSE}"
echo "[INFO] Log directory: ${LOG_DIR}"

# Define the suffix of the log file to be merged.（Order matters）
LOG_TYPES=(
  "trim"
  "filter"
  "count_unique"
  "map_genome"
  "map_rRNA_5S"
  "map_ncRNA_nuclear_tRNA"
  "process_rna_seq"
  "calculate_rpm"
)

# Label
LABEL_TYPES=(
  "trim"
  "filter"
  "count_unique"
  "map_genome"
  "map_rRNA"
  "map_ncRNA"
  "process_rna_seq"
  "calculate_rpm"
)

##### Merge {sample} logs
for SAMPLE in ${SAMPLES}; do
  echo "[INFO] Merging logs for sample: ${SAMPLE}"

  OUT_LOG="${LOG_DIR}/${SAMPLE}_rsRNAQuant.log"
  : > "${OUT_LOG}"  


  for i in "${!LOG_TYPES[@]}"; do
    LOG_TYPE="${LOG_TYPES[$i]}"
    LABEL="${LABEL_TYPES[$i]}"

    LOG_FILE="${LOG_DIR}/${SAMPLE}_${LOG_TYPE}.log"

    if [[ "$NCRNA_FLAG" == "false" && "$LOG_TYPE" == "map_ncRNA" ]]; then
        # Skip the log of ncRNA
        continue
    else
      echo "===============================" >> "${OUT_LOG}"
      echo "[${SAMPLE}] ${LABEL}"           >> "${OUT_LOG}"
      echo "===============================" >> "${OUT_LOG}"

      if [[ -f "${LOG_FILE}" ]]; then
        cat "${LOG_FILE}" >> "${OUT_LOG}"
      else
        echo "[WARNING] Missing log file: ${LOG_FILE}" >> "${OUT_LOG}"
      fi

    fi

    echo -e "\n" >> "${OUT_LOG}"
  done
done

# Delete
for SAMPLE in ${SAMPLES}; do
  for LOG_TYPE in "${LOG_TYPES[@]}"; do
    rm -f "${LOG_DIR}/${SAMPLE}_${LOG_TYPE}.log"
  done
  rm -f "${LOG_DIR}/${SAMPLE}_map_ncRNA_"*.log
  rm -f "${LOG_DIR}/${SAMPLE}_map_rRNA_"*.log
done







##### Merge {GSE} logs
GSE_LOG_TYPES=(
  "align_summary"
  "calculate_CountMatrix"
  "calculate_ExprMatrix"
  "plot_length"
  "plot_coverage"
  "plot_DA"
)

##### Merge {sample} logs
echo "[INFO] Merging logs for project: ${GSE}"

OUT_LOG="${LOG_DIR}/${GSE}_rsRNAQuant.log"
: > "${OUT_LOG}"   

for GSE_LOG_TYPE in "${GSE_LOG_TYPES[@]}"; do
  LOG_FILE="${LOG_DIR}/${GSE}_${GSE_LOG_TYPE}.log"

  echo "===============================" >> "${OUT_LOG}"
  echo "[${GSE}] ${GSE_LOG_TYPE}"        >> "${OUT_LOG}"
  echo "===============================" >> "${OUT_LOG}"

  if [[ -f "${LOG_FILE}" ]]; then
    cat "${LOG_FILE}" >> "${OUT_LOG}"
  else
    echo "[WARNING] Missing log file: ${LOG_FILE}" >> "${OUT_LOG}"
  fi

  echo -e "\n" >> "${OUT_LOG}"
done

##### Delete 
for GSE_LOG_TYPE in "${GSE_LOG_TYPES[@]}"; do
  rm -f "${LOG_DIR}/${GSE}_${GSE_LOG_TYPE}.log"
done


# echo "Cleaning temporary directories..." >&2
rm -f -r "${RESULTS_DIR}/${GSE}/filtered/"
rm -f -r "${RESULTS_DIR}/${GSE}/trimmed/"
rm -f "${RESULTS_DIR}/${GSE}/plots/${GSE}*.pdf"

for SAMPLE in ${SAMPLES}; do
  rm -f "${RESULTS_DIR}/${GSE}/stats/${SAMPLE}_rRNA_output.txt"
done

# Only keep rsRNA files
if [ "$NCRNA_FLAG" == "false" ]; then
  for SAMPLE in ${SAMPLES}; do
    rm -f "${RESULTS_DIR}/${GSE}/stats/${SAMPLE}_summary.txt"
    rm -f "${RESULTS_DIR}/${GSE}/stats/${SAMPLE}_output.txt"
    rm -f "${RESULTS_DIR}/${GSE}/stats/${SAMPLE}_length_distribution.txt"
  done
fi

<<EOF
if [ "$FILES_FLAG" = "remove" ]; then
  rm "${RESULTS_DIR}/${GSE}/plots/${GSE}*_Data.txt" 
  rm "${RESULTS_DIR}/${GSE}/plots/${GSE}_rRNA_End_RPM.txt"
  rm "${RESULTS_DIR}/${GSE}/plots/${GSE}_rRNA_subtype_RPM.txt"
  rm "${RESULTS_DIR}/${GSE}/plots/${GSE}_rRNA_Length_distribution_RPM.txt"
fi
EOF

if [ "$FILES_FLAG" = "remove" ]; then
  rm -f -r "${RESULTS_DIR}/${GSE}/mapped/"
fi


##### Generate a flag file for Snakemake to track.
merge_flag="${LOG_DIR}/.all.done"
touch "${merge_flag}"

echo "[INFO] Log merge finished. Flag created: ${merge_flag}"
