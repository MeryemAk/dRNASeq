#!/bin/bash

INDIR="$HOME/Shared/00-stage/dRNASeq2/1.data/merged"
OUTDIR="$HOME/Shared/00-stage/dRNASeq2/3.trimming/porechop"
THREADS=4

mkdir -p "$OUTDIR"

for file in "$INDIR"/*.fastq; do
    echo "Trimming with Porechop for $file..."

    BASENAME=$(basename "$file" .fastq)
    FILE_NAME_ONLY=$(basename "$file")

    porechop -t ${THREADS} -i "$file" -o "$OUTDIR/${BASENAME}_trimmed.fastq"


    if [ $? -ne 0 ]; then
        echo "Trimming failed for ${FILE_NAME_ONLY}"
        exit 1
    fi
    echo "Trimming completed for ${FILE_NAME_ONLY}"
done
echo "Done trimming with Porechop"