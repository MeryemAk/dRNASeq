#!/bin/bash

# This master script orchestrates the execution of all dRNASeq processing scripts.
# It assumes all individual scripts and the 'config.conf' file are in the same directory as this script.

# Exit immediately if a command exits with a non-zero status.
set -e

echo "Starting the dRNASeq processing pipeline..."
echo "--------------------------------------------------"

# --- Step 1: Linux Setup and Conda Environment Installation ---
# This script should ideally be run once to set up the environment.
# It will exit if the conda environment already exists.
echo "Running 1.linuxsetup.sh to set up the Conda environment..."
bash 1.linuxsetup.sh
echo "1.linuxsetup.sh completed."
echo "--------------------------------------------------"

# After linuxsetup.sh, you might need to restart your terminal or
# manually activate the conda environment if it wasn't already active.
echo "Activating dRNAseq conda environment..."
# Check if conda is initialized and source it if not already
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
else
    echo "Error: Conda initialization script not found. Please ensure Miniconda is installed correctly."
    exit 1
fi

# Activate the dRNAseq conda environment
conda activate dRNAseq
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate 'dRNAseq' conda environment. Please check your setup."
    exit 1
fi
echo "dRNAseq conda environment activated."
echo "--------------------------------------------------"

# --- Step 2: Merge FASTQ files ---
echo "Running 2.merge.sh to merge FASTQ files..."
bash 2.merge.sh
echo "2.merge.sh completed."
echo "--------------------------------------------------"

# --- Step 3: Filter rRNA ---
echo "Running 3.filter_rRNA.sh to filter rRNA..."
bash 3.filter_rRNA.sh
echo "3.filter_rRNA.sh completed."
echo "--------------------------------------------------"

# --- Step 4: Trim FASTQ files with Pychopper ---
echo "Running 4.pychopper.sh to trim FASTQ files..."
bash 4.pychopper.sh
echo "4.pychopper.sh completed."
echo "--------------------------------------------------"

# --- Step 5: Perform QC with NanoPack ---
echo "Running 5.qc.sh to perform QC analysis..."
bash 5.qc.sh
echo "5.qc.sh completed."
echo "--------------------------------------------------"

# --- Step 6: Map reads with Minimap2 ---
echo "Running 6.mapping.sh to map reads..."
bash 6.mapping.sh
echo "6.mapping.sh completed."
echo "--------------------------------------------------"

# --- Step 7: Count reads with Bambu Runner ---
echo "Running 7.counting.sh to count reads..."
bash 7.counting.sh
echo "7.counting.sh completed."
echo "--------------------------------------------------"

# --- Step 8: Kraken2 Classification ---
echo "Running 8.kraken.sh for Kraken2 classification..."
bash 8.kraken.sh
echo "8.kraken.sh completed."
echo "--------------------------------------------------"

echo "All dRNASeq processing pipeline steps completed successfully!"