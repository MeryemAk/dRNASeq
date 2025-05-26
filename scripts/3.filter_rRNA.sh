#!/bin/bash

# Source the configuration file
source config.conf

# --- Handle Input ---
echo "Input directory for filtering: $FILTER_RRNA_INPUT_DIR"
echo "Output directory for filtered files: $FILTER_RRNA_OUTPUT_DIR"

# Create output directory if it doesn't exist
mkdir -p "$FILTER_RRNA_OUTPUT_DIR"

# Loop through all FASTQ files in the input directory
for FILE in "$FILTER_RRNA_INPUT_DIR"/*.fastq; do
    BASENAME=$(basename "$FILE" .fastq)
    SAM_FILE="$FILTER_RRNA_OUTPUT_DIR/${BASENAME}_mapped.sam"
    FILTERED_FASTQ="$FILTER_RRNA_OUTPUT_DIR/${BASENAME}_filtered.fastq"

    echo "Processing $FILE..."
    minimap2 -ax map-ont "$FILTER_RRNA_DB" "$FILE" > "$SAM_FILE"

    echo "Extracting unmapped reads..."
    samtools fastq -f 4 "$SAM_FILE" > "$FILTERED_FASTQ"

    # Delete the SAM file after extracting unmapped reads
    rm -f "$SAM_FILE"
    
done

echo "Filtering completed!"
