#!/usr/bin/env python3
# This script takes an accession file containing SRR IDs and downloads the corresponding zipped FASTQ files.
# The script supports both single-end (SE) and paired-end (PE) reads.

import os
import sys # for command line arguments

#########################################################################################################
# USAGE
#########################################################################################################
def usage():
    print(
        "Usage: python3 2.1_SRR_download.py <accession_file> <SE|PE> [work_dir]\n"
        "\n"
        "Parameters:\n"
        "  <accession_file> : Path to the file containing SRR accession numbers (one per line).\n"
        "  <SE|PE>          : Specify 'SE' for single-end reads or 'PE' for paired-end reads.\n"
        "  [work_dir]       : Path to the output directory.\n"
    )
    sys.exit(1)
#########################################################################################################
# INPUT VALIDATION
#########################################################################################################
# Check if the correct number of arguments is provided
if len(sys.argv) < 3:
    usage()

# Extract information from provided arguments
accession_file = sys.argv[1]
read_type = sys.argv[2].upper()  # Convert to uppercase for consistency
work_dir = sys.argv[3]

# Validate the read type
if read_type not in ["SE", "PE"]:
    print("ERROR: Invalid read type. Use 'SE' for single-end or 'PE' for paired-end.")
    sys.exit(1)

# Check if the accession file exists
if not os.path.exists(accession_file):
    print(f"ERROR: Accession file '{accession_file}' not found.")
    sys.exit(1)

# Check if output directory exists, if not create it
if not os.path.exists(work_dir):
    os.makedirs(work_dir)

#########################################################################################################
# READ ACCESSION FILE
#########################################################################################################
# Create a dictionary of samples from the accession file
samples = {}
with open(accession_file, "r") as f:
    for line in f:
        line = line.strip()
        if line:
            sample_id = f"Sample_{len(samples) + 1}"  # Generate a sample ID
            samples[sample_id] = line

#########################################################################################################
# RUN FASTQ-DUMP
#########################################################################################################
# downloading each given accession number
for sample_id, accession in samples.items():
    print(f'--- Currently downloading:{accession} ---')

    # Construct the fastq-dump command based on the read type
    if read_type == "PE":
        cmd_fastqdump = f'fastq-dump --gzip --split-3 -O {work_dir} {accession}'
    else:  # SE
        cmd_fastqdump = f'fastq-dump --gzip -O {work_dir} {accession}'

    # Show the command
    print(f"Running command: {cmd_fastqdump}")

    # Execute the command
    exit_code = os.system(cmd_fastqdump)

    # Check if the command was successful
    if exit_code != 0:
        print(f"ERROR: fastq-dump failed for accession {accession}.")
    else:
        print(f"Successfully downloaded {accession}.\n")