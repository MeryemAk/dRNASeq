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
mkdir -p "${KRAKEN_OUTPUT_DIR}"

# Check if Kraken2 database exists
if [ ! -d "$KRAKEN_DB" ]; then
    echo "Error: Kraken2 database not found at $KRAKEN_DB"
    exit 1
fi

echo "-------------------------"

# Loop through all sample subdirectories in the mapping output directory
find "$KRAKEN_INPUT_DIR" -maxdepth 1 -type d -not -path "$KRAKEN_INPUT_DIR" | while read SAMPLE_DIR; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_filtered_trimmed_bacteria_unmapped.fastq"

    # Check if the unmapped FASTQ file exists
    if [ ! -f "$FILE" ]; then
        echo "Warning: Unmapped FASTQ file not found: $FILE"
        continue
    fi

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