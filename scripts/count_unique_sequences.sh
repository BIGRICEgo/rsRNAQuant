#!/bin/bash
# Usage: bash count_unique_sequences.sh <filtered_fastq> <unique_fa> <len_dist> <tmp_seq> <tmp_count> <log>

FILTERED_FASTQ="$1"
UNIQUE_FA="$2"
LEN_DIST="$3"
TMP_SEQ="$4"
TMP_COUNT="$5"
LOG="$6"


echo "Starting counting unique sequences." >> "$LOG"
awk 'NR%4==2' "$FILTERED_FASTQ" > "$TMP_SEQ"
sort "$TMP_SEQ" | uniq -c | sort -k1,1nr > "$TMP_COUNT"
awk 'BEGIN{OFS=""}{printf(">t%08d\t%d\n%s\n", NR, $1, $2)}' "$TMP_COUNT" > "$UNIQUE_FA"
awk '{print length($2), $1}' "$TMP_COUNT" | \
    awk '{a[$1]+=$2}END{for(i in a) print i, a[i]}' | sort -n > "$LEN_DIST"


rm -f "$TMP_SEQ" "$TMP_COUNT"
echo "Successfully completed counting unique sequences." >> "$LOG"

