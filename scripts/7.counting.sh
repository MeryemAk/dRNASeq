#!/bin/bash

# Source the configuration file
source config.conf

# --- Handle Input ---
echo "Counting with Bambu Runner Docker..."
echo "Output directory: $COUNTING_OUTPUT_DIR"
echo "Processing data from: $COUNTING_INPUT_DIR"
echo "Species list: $COUNTING_SPECIES_LIST"
echo "Annotation files: $COUNTING_ANNOTATIONS"
echo "Genome files: $COUNTING_GENOMES"

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
        GENOME_FILE="${COUNTING_GENOMES[$i]}"
        
        echo "Processing species: $SPECIES using annotation: $ANNOTATION_FILE"

        # Find the BAM file for the species in the sample directory
        BAM_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_mapped_${SPECIES}_sorted.bam"
        if [ ! -f "$BAM_FILE" ]; then
            echo "Warning: BAM file not found for species $SPECIES in sample $SAMPLE_NAME. Skipping..."
            continue
        fi

        echo "Running Bambu Runner for sample: $SAMPLE_NAME, species: $SPECIES"

        # Run the Bambu Runner Docker container
        docker run --rm \
            -v "${SAMPLE_DIR}:/data" \
            -v "${COUNTING_OUTPUT_DIR}:/output" \
            -v "$HOME/dRNASeq/reference_genomes:/reference_genomes" \
            mathiasverbeke/bambu_runner:latest \
            run_bambu.R \
            --reads "/data/$(basename "$BAM_FILE")" \
            --annotations "/reference_genomes/$(basename "$ANNOTATION_FILE")" \
            --genome "/reference_genomes/$(basename "$GENOME_FILE")" \
            --output-dir "/output/${SAMPLE_NAME}" \
            --ncore 1 \
            --stranded no \
            --quant yes \
            --discovery no \
            --verbose yes \
            --rc-out-dir "/output/${SAMPLE_NAME}/rc_cache" \
            --low-memory no

        echo "Finished processing: $SAMPLE_NAME - $SPECIES"
    done
done

echo "All samples quantification completed!"
