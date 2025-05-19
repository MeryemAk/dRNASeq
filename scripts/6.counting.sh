#!/bin/bash

# Source the configuration file
source config.conf

# --- Handle Input ---
echo "Counting with NanoCount..."
echo "Output directory: $COUNTING_OUTPUT_DIR"
echo "Processing data from: $COUNTING_INPUT_DIR"
echo "Species list: $COUNTING_SPECIES_LIST"
echo "Annotation fils: $COUNTING_ANNOTATIONS"

# Create the output directory if it doesn't exist
mkdir -p "${COUNTING_OUTPUT_DIR}"

echo "-------------------------"

# Loop through each sample directory in the mapping output directory
for SAMPLE_DIR in "$COUNTING_INPUT_DIR"/*; do
    # Skip if not a directory
    [ -d "$SAMPLE_DIR" ] || continue

    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    SAMPLE_OUTPUT_DIR="${COUNTING_OUTPUT_DIR}/${SAMPLE_NAME}"
    mkdir -p "$SAMPLE_OUTPUT_DIR"

    echo "Processing sample: $SAMPLE_NAME"

    # Loop through each species by index
    for i in "${!COUNTING_SPECIES_LIST[@]}"; do
        SPECIES="${COUNTING_SPECIES_LIST[$i]}"
        ANNOTATION_FILE="${COUNTING_ANNOTATIONS[$i]}"
        
        echo "Processing species: $SPECIES using annotation: $ANNOTATION_FILE"

        # Find the BAM file for the species in the sample directory
        BAM_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_mapped_${SPECIES}_sorted.bam"
        if [ ! -f "$BAM_FILE" ]; then
            echo "Warning: BAM file not found for species $SPECIES in sample $SAMPLE_NAME. Skipping..."
            continue
        fi

        OUT_FILE="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME}_${SPECIES}_nanocount.tsv"
        
        echo "Quantifying sample: $SAMPLE_NAME for species: $SPECIES"
        
        # Run NanoCount using the BAM file and corresponding annotation (GTF/GFF)
        NanoCount -i "$BAM_FILE" -a "$ANNOTATION_FILE" -o "$OUT_FILE"
        
        echo "Finished processing: $SAMPLE_NAME - $SPECIES"
    done
done

echo "All samples quantification completed!"
