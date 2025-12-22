#!/bin/bash
# align_summary.sh
# Usage: bash align_summary.sh <GSE> <results_dir> <samples> <NCRNA_FLAG>

GSE="$1"
RESULTS_DIR="$2"
SAMPLES="$3"  # Space-separated sample names
NCRNA_FLAG="$4" 

rRNA_type_order=('5S' '5.8S' '12S' '16S' '18S' '28S' '45S')
ncRNA_type_order=('nuclear_tRNA' 'mt_RNA' 'yRNA' 'miRNA' 'snoRNA' 'snRNA' 'vaultRNA' 'piRNA')

output_dir="${RESULTS_DIR}/stats"
mkdir -p "$output_dir"
output_file="${output_dir}/${GSE}_alignment_summary.csv"

# Header
echo -n "SAMPLE,Raw,Preprocess,Filter,Genome," > "$output_file"
for rRNA in "${rRNA_type_order[@]}"; do echo -n "${rRNA}_aligned," >> "$output_file"; done
if [ "$NCRNA_FLAG" = "true" ]; then
    for ncRNA in "${ncRNA_type_order[@]}"; do echo -n "${ncRNA}_aligned," >> "$output_file"; done
fi
echo >> "$output_file"


for SAMPLE in $SAMPLES; do
    echo -n "${SAMPLE}," >> "$output_file"

    # Raw
    stat_file="${RESULTS_DIR}/trimmed/${SAMPLE}.fastq.gz_trimming_report.txt"
    if [ -f "$stat_file" ]; then
        reads=$(grep "Total reads processed:" "$stat_file" | awk -F':' '{print $2}' | tr -d '[:space:],')
        echo -n "$reads," >> "$output_file"
    else
        echo -n "0," >> "$output_file"
    fi

    # Preprocess
    stat_file="${RESULTS_DIR}/filtered/${SAMPLE}_match_univec.bowtie.stat"
    if [ -f "$stat_file" ]; then
        reads=$(grep "reads processed:" "$stat_file" | awk -F':' '{print $2}')
        echo -n "$reads," >> "$output_file"
    else
        echo -n "0," >> "$output_file"
    fi

    # Filter
    if [ -f "$stat_file" ]; then
        reads=$(grep "reads that failed to align:" "$stat_file" | awk -F':' '{print $2}')
        echo -n "$reads," >> "$output_file"
    else
        echo -n "0," >> "$output_file"
    fi

    # Genome
    stat_file="${RESULTS_DIR}/mapped/${SAMPLE}_match_genome.bowtie.stat"
    if [ -f "$stat_file" ]; then
        reads=$(grep "reads with at least one alignment:" "$stat_file" | awk -F':' '{print $2}')
        echo -n "$reads," >> "$output_file"
    else
        echo -n "0," >> "$output_file"
    fi

    # rRNA
    for rRNA in "${rRNA_type_order[@]}"; do
        stat_file="${RESULTS_DIR}/mapped/${SAMPLE}_MG_match_${rRNA}_rRNA.bowtie.stat"
        if [ -f "$stat_file" ]; then
            reads=$(grep "reads with at least one alignment:" "$stat_file" | awk -F':' '{print $2}')
            echo -n "$reads," >> "$output_file"
        else
            echo -n "0," >> "$output_file"
        fi
    done

if [ "$NCRNA_FLAG" = "true" ]; then
    # ncRNA
    for ncRNA in "${ncRNA_type_order[@]}"; do
        stat_file="${RESULTS_DIR}/mapped/${SAMPLE}_MG_match_${ncRNA}_ncRNA.bowtie.stat"
        if [ -f "$stat_file" ]; then
            reads=$(grep "reads with at least one alignment:" "$stat_file" | awk -F':' '{print $2}')
            echo -n "$reads," >> "$output_file"
        else
            echo -n "0," >> "$output_file"
        fi
    done
fi

    echo >> "$output_file"
done

echo "Alignment summary written to ${output_file}"

