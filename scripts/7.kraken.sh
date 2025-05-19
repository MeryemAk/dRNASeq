#!/bin/bash

# Source the configuration file
source config.conf

# --- Handle Input ---
echo "Kraken2 classification..."
echo "Output directory: $KRAKEN_OUTPUT_DIR"
echo "Processing data from: $KRAKEN_INPUT_DIR"
echo "Kraken2 database: $KRAKEN_DB"
echo "Number of threads: $KRAKEN_THREADS"

# Create the output directory if it doesn't exist
mkdir -p "${COUNTING_OUTPUT_DIR}"

# Check if Kraken2 database exists
if [ ! -d "$KRAKEN_DB" ]; then
    echo "Error: Kraken2 database not found at $KRAKEN_DB"
    exit 1
fi

echo "-------------------------"

# Loop through all FASTQ files in the input directory
for FILE in "$KRAKEN_INPUT_DIR"/*.fastq; do
    SAMPLE_NAME=$(basename "$FILE" _merged.fastq)
    SAMPLE_OUTPUT_DIR="${KRAKEN_OUTPUT_DIR}/${SAMPLE_NAME}"

    # Create a subdirectory for each sample
    mkdir -p "$SAMPLE_OUTPUT_DIR"

    echo "Starting Kraken2 classification for $SAMPLE_NAME"

    # Run Kraken2
    kraken2 --db "$KRAKEN_DB" \
        --threads "$KRAKEN_THREADS" \
        --report "${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME}_report.txt" \
        --output "${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME}_output.txt" \
        --use-names \
        --classified-out "${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME}_classified.fastq" \
        --unclassified-out "${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME}_unclassified.fastq" \
        "$FILE"

    echo "Finished processing sample: $SAMPLE_NAME"
    echo "-------------------------"
done

echo "All samples processed successfully!"