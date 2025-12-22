#!/bin/bash
# Usage: bash map_rRNA.sh <sample> <results_dir> <annotations_dir> <threads>

SAMPLE="$1"
RESULTS_DIR="$2"
ANNOT_DIR="$3"
THREADS="$4"

echo "Running rRNA alignment for ${SAMPLE}..."

input_file="${RESULTS_DIR}/mapped/${SAMPLE}_genome_mapped.fa"
rRNA_types=('5S' '5.8S' '12S' '16S' '18S' '28S' '45S')

for i in "${!rRNA_types[@]}"; do
    rRNA=${rRNA_types[$i]}

    if [ $i -eq 0 ]; then
        in_file="${input_file}"
    else
        prev=${rRNA_types[$((i-1))]}
        in_file="${RESULTS_DIR}/mapped/${SAMPLE}_MG_unmatch_${prev}_rRNA.fa"
    fi

    out_match="${RESULTS_DIR}/mapped/${SAMPLE}_MG_match_${rRNA}_rRNA.fa"
    out_unmatch="${RESULTS_DIR}/mapped/${SAMPLE}_MG_unmatch_${rRNA}_rRNA.fa"
    out_detail="${RESULTS_DIR}/mapped/${SAMPLE}_detail_match_${rRNA}_rRNA"
    out_stat="${RESULTS_DIR}/mapped/${SAMPLE}_MG_match_${rRNA}_rRNA.bowtie.stat"
    index="${ANNOT_DIR}/rRNAdb/human_rRNA_${rRNA}"

    echo "  Mapping to ${rRNA} ..."
    bowtie -f "${in_file}" -v 1 -a -p "${THREADS}" --fullref --norc \
           --al "${out_match}" --un "${out_unmatch}" \
           -x "${index}" > "${out_detail}" 2> "${out_stat}"
done

echo "Finished rRNA mapping for ${SAMPLE}"
