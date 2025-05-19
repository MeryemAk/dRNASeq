##!/bin/bash
# Script to process FASTQ files using NanoPack for QC
# Results will be placed directly in the 2.qc folder.

# Source the configuration file
source config.conf

# --- Handle Input ---
echo "Starting QC analysis on merged data using NanoPack..."
echo "Output directory: $QC_OUTPUT_DIR"
echo "Processing data from: $MERGED_DATA_DIR"
echo "Number of threads (QC): $QC_THREADS"

# Create the output directory if it doesn't exist
mkdir -p "${QC_OUTPUT_DIR}"

echo "-------------------------"

# --- Process FASTQ Files ---
# Find all FASTQ files within the merged data directory
find "$MERGED_DATA_DIR" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | while IFS= read -r fq_path; do
  echo "Processing: $fq_path"

  # --- NanoPlot QC using NanoPack ---
  fq_base=$(basename "$fq_path")
  fq_name="${fq_base%.*}" # Remove extension

  echo "    Running NanoPlot..."
  NanoPlot -t "$DEFAULT_THREADS_QC" -o "${QC_OUTPUT_DIR}" --prefix "${fq_name}_" --N50 --fastq "$fq_path"

  echo "  Finished processing: $fq_path"
  echo "-------------------------"
done

# --- Run NanoComp after all NanoPlot executions ---
echo "Running NanoComp on all processed FASTQ files..."
NanoComp -t "$DEFAULT_THREADS_QC" -o "${QC_OUTPUT_DIR}" -f png --plot violin --fastq "$MERGED_DATA_DIR"/*.fastq

echo "QC analysis using NanoPack complete. Results are in: $QC_OUTPUT_DIR"

# --- Clean up ---

