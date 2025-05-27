#!/bin/bash
# Script to trim FASTQ files using Pychopper

# Source the configuration file
source config.conf

# --- Handle Input ---
echo "Trimming FASTQ files with Pychopper..."
echo "Output directory: $TRIMMING_OUTPUT_DIR"
echo "Processing data from: $TRIMMING_INPUT_DIR"
echo "Number of threads: $TRIMMING_THREADS"

# Create the output directory if it doesn't exist
mkdir -p "${TRIMMING_OUTPUT_DIR}"

echo "-------------------------"

mkdir -p "$TRIMMING_OUTPUT_DIR"

for file in "$TRIMMING_INPUT_DIR"/*.fastq; do
    echo "Trimming with Pychopper for $file..."

    BASENAME=$(basename "$file" .fastq)
    FILE_NAME_ONLY=$(basename "$file")

    pychopper -t $TRIMMING_THREADS "$file" "$TRIMMING_OUTPUT_DIR/${BASENAME}_trimmed.fastq" -r "pychopper_report.pdf"

    if [ $? -ne 0 ]; then
        echo "Trimming failed for ${FILE_NAME_ONLY}"
        exit 1
    fi
    echo "Trimming completed for ${FILE_NAME_ONLY}"

    echo "-------------------------"
done
echo "Done trimming with Pychopper"