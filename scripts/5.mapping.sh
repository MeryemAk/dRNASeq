#!/bin/bash

# Source the configuration file
source config.conf

# --- Handle Input ---
echo "Mapping with Minimap2..."
echo "Output directory: $MAPPING_OUTPUT_DIR"
echo "Processing data from: $MAPPING_INPUT_DIR"
echo "Species list: $MAPPING_SPECIES_LIST"
echo "Reference indexes: $MAPPING_REF_INDEXES"
echo "Mapping parameters: $MAPPING_PARAMS"

# Create the output directory if it doesn't exist
mkdir -p "${MAPPING_OUTPUT_DIR}"


# Loop over FASTQ files in MAPPING_INPUT_DIR
for input_file in "$MAPPING_INPUT_DIR"/*; do
    # Only process files (skip directories, if any)
    [ -f "$input_file" ] || continue

    # Get the sample name (remove the file extension)
    sample_name=$(basename "$input_file")
    sample_name="${sample_name%.*}"
    echo "===== Processing sample: $sample_name ====="
    
    # Start with the original FASTQ file as the input for the first mapping
    previous_unmapped="$input_file"
    
    # Create the output directory for the sample if it doesn't exist
    sample_outdir="${MAPPING_OUTPUT_DIR}/${sample_name}"
    if [ ! -d "$sample_outdir" ]; then
        mkdir -p "$sample_outdir"
        echo "Created directory: $sample_outdir"
    fi

    # Loop through each species
    for i in "${!MAPPING_SPECIES_LIST[@]}"; do
        species="${MAPPING_SPECIES_LIST[$i]}"
        ref_index="${MAPPING_REF_INDEXES[$i]}"
        params="${MAPPING_PARAMS[$i]}"

        echo "=== Mapping to ${species^} reference genome ==="
        
        # Check if the reference genome index exists
        if [ ! -f "$ref_index" ]; then
            echo "Error: ${species^} reference genome file not found: $ref_index"
            exit 1
        else
            echo "${species^} reference genome file found: $ref_index"
        fi

        # Map reads with minimap2
        sam_output="${sample_outdir}/${sample_name}_${species}.sam"
        cmd="minimap2 ${params} ${ref_index} ${previous_unmapped} > ${sam_output}"
        echo "Mapping command: $cmd"
        eval "$cmd"
        echo "Mapped ${previous_unmapped} to ${species^} genome: ${sam_output}"

        # Extract unmapped reads using samtools
        echo "===== Extracting unmapped reads from ${species^} mapping ====="
        unmapped_output="${sample_outdir}/${sample_name}_${species}_unmapped.fastq"
        cmd_unmapped="samtools fastq -f 4 ${sam_output} > ${unmapped_output}"
        echo "Extraction command: $cmd_unmapped"
        eval "$cmd_unmapped"
        echo "Unmapped reads saved as: ${unmapped_output}"

        # Sort and index the mapped reads using samtools
        echo "===== Sorting and indexing mapped reads for ${species^} ====="
        sorted_bam="${sample_outdir}/${sample_name}_mapped_${species}_sorted.bam"
        cmd_sort="samtools sort ${sam_output} -o ${sorted_bam}"
        echo "Sort command: $cmd_sort"
        eval "$cmd_sort"
        
        cmd_index="samtools index ${sorted_bam}"
        echo "Index command: $cmd_index"
        eval "$cmd_index"
        echo "Sorted and indexed BAM: ${sorted_bam}"

        # Generate alignment statistics using seqkit
        echo "===== Generating alignment statistics for ${species^} ====="
        stats_output="${sample_outdir}/${sample_name}_${species}_stats.txt"
        cmd_stats="seqkit bam -s ${sorted_bam} > ${stats_output}"
        echo "Alignment stats command: $cmd_stats"
        eval "$cmd_stats"
        echo "Alignment statistics saved as: ${stats_output}"
        
        # Update the input for the next iteration to be the unmapped reads from current mapping
        previous_unmapped="${unmapped_output}"
        echo "--------------------------------------------------"
    done
    echo "===================================================="
done

echo "Mapping completed for all samples."
