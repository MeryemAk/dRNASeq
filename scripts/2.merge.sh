#!/bin/bash
# Script to merge fastq files based on directory structure.

# Source the configuration file
source config.conf

# --- Handle Input ---
echo "Input directory for merging: $MERGE_INPUT_DIR"
echo "Output directory for merged files: $MERGE_OUTPUT_DIR"

# Create the output directory if it doesn't exist
mkdir -p "${MERGE_OUTPUT_DIR}"

echo "-------------------------"

# Find all barcode directories (e.g., barcode01, barcode02, ...)
find "${MERGE_INPUT_DIR}" -maxdepth 1 -type d -name "barcode*" | while IFS= read -r barcode_dir; do
  # Extract the barcode name from the directory path.
  barcode_name=$(basename "${barcode_dir}")

  # Define the output file name for the merged files of this barcode.
  output_file="${MERGE_OUTPUT_DIR}/${barcode_name}_merged.fastq"

  # Check if the output file already exists
  if [ -e "${output_file}" ]; then
    echo "Output file ${output_file} already exists. Overwriting."
  fi

  # Use cat to concatenate all FASTQ files in the barcode directory.
  find "${barcode_dir}" -type f -name "*.fastq" -print0 | sort -z | xargs -0 cat > "${output_file}"

  if [ $? -eq 0 ]; then
    echo "Merged all FASTQ files in ${barcode_dir} into ${output_file}"
  else
    echo "Error merging FASTQ files in ${barcode_dir}."
  fi

  echo "-------------------------"
done

echo "Finished merging FASTQ files."