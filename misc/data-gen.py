#!/usr/bin/env python3

import os
import random

# Input FASTA files
input_files = ["SRR26937771.fastq", "SRR23886070.fastq", "SRR29716937.fastq", "SRR20324764.fastq"]

# Max number of reads to sample per file
# This defines the maximum number of reads to randomly sample from each input file
max_reads_per_file = 100

# Output directory
output_dir = "data"
os.makedirs(output_dir, exist_ok=True) # If output directory doesn't exist, create it

# Generate 10 output files
for i in range(1, 11): # Generate 10 output files named combined_1.fastq to combined_10.fastq
    combined_records = []

    # Process each input file
    for file in input_files:
        with open(file, 'r') as f:
            lines = f.readlines()

        # Group lines into FASTQ records (4 line format)
        # Each FASTQ record consists of 4 lines: header, sequence, separator, and quality
        records = [lines[j:j+4] for j in range(0, len(lines), 4)]
        total_available = len(records)

        # Determine the number of records to sample from this file
        # Randomly choose a number between 1 and the smaller of max_reads_per_file or total_available
        n_to_sample = random.randint(1, min(max_reads_per_file, total_available))

        # Randomly sample the determined number of records
        sampled = random.sample(records, n_to_sample)

        # Append the sampled records to the combined list
        combined_records.extend(sampled)

    # Shuffle the combined records to randomize their order
    # This ensures that records from different files are mixed together
    random.shuffle(combined_records)

    # Flatten the list of records into a single list of lines
    # Each record (4 lines) is unpacked into individual lines
    output_lines = [line for record in combined_records for line in record]

    # Define the output file path
    output_path = os.path.join(output_dir, f"combined_{i}.fastq")

    # Write the combined and shuffled records to the output file
    with open(output_path, 'w') as out_f:
        out_f.writelines(output_lines)

    # Print a message indicating the output file has been generated
    print(f"Generated {output_path} with {len(combined_records)} sequences")
