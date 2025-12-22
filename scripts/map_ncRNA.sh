#!/bin/bash
# Usage: bash map_ncRNA.sh <sample> <results_dir> <annotations_dir> <threads>

SAMPLE="$1"
RESULTS_DIR="$2"
ANNOT_DIR="$3"
THREADS="$4"

echo "Running ncRNA alignment for ${SAMPLE}..."

input_file="${RESULTS_DIR}/mapped/${SAMPLE}_MG_unmatch_45S_rRNA.fa"

ncRNA_types=('nuclear_tRNA' 'mt_RNA' 'yRNA' 'miRNA' 'snoRNA' 'snRNA' 'vaultRNA' 'piRNA')
ncRNA_indices=(
    "${ANNOT_DIR}/GtRNAdb/hg38-tRNAs"
    "${ANNOT_DIR}/GtRNAdb/hg38-mt_tRNAs"
    "${ANNOT_DIR}/YRNA/human_YRNA"
    "${ANNOT_DIR}/miRBase/miRBase_22-hsa"
    "${ANNOT_DIR}/snoDB/snoDB"
    "${ANNOT_DIR}/RNAcentral/human_snRNA_Refseq"
    "${ANNOT_DIR}/RNAcentral/vaultRNA"
    "${ANNOT_DIR}/piRBase/piRBase_hsa.gold"
)

for i in "${!ncRNA_types[@]}"; do
    ncRNA=${ncRNA_types[$i]}
    index=${ncRNA_indices[$i]}

    if [ $i -eq 0 ]; then
        in_file="${input_file}"
    else
        prev=${ncRNA_types[$((i-1))]}
        in_file="${RESULTS_DIR}/mapped/${SAMPLE}_MG_unmatch_${prev}_ncRNA.fa"
    fi

    out_match="${RESULTS_DIR}/mapped/${SAMPLE}_MG_match_${ncRNA}_ncRNA.fa"
    out_unmatch="${RESULTS_DIR}/mapped/${SAMPLE}_MG_unmatch_${ncRNA}_ncRNA.fa"
    out_detail="${RESULTS_DIR}/mapped/${SAMPLE}_detail_match_${ncRNA}_ncRNA"
    out_stat="${RESULTS_DIR}/mapped/${SAMPLE}_MG_match_${ncRNA}_ncRNA.bowtie.stat"

    echo "  Mapping to ${ncRNA} ..."
    bowtie -f "${in_file}" -v 1 -k 3 -p "${THREADS}" --fullref --norc --best --strata \
           --al "${out_match}" --un "${out_unmatch}" \
           -x "${index}" > "${out_detail}" 2> "${out_stat}"
done

echo "Finished ncRNA mapping for ${SAMPLE}"
