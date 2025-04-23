#!/usr/bin/env python3

import os
import random

# Input FASTA files
input_files = ["SRR26937771.fastq", "SRR23886070.fastq", "SRR29716937.fastq", "SRR20324764.fastq"]

# Max number of reads to sample per file
max_reads_per_file = 100

# Output directory
output_dir = "combined_fasta"
os.makedirs(output_dir, exist_ok=True)

for i in range(1, 11):  # Generate 10 output files
    combined_records = []

    for file in input_files:
        with open(file, 'r') as f:
            lines = f.readlines()

        # Group lines into FASTQ records (4 line format)
        records = [lines[j:j+2] for j in range(0, len(lines), 2)]
        total_available = len(records)
        n_to_sample = random.randint(1, min(max_reads_per_file, total_available))
        sampled = random.sample(records, n_to_sample)
        combined_records.extend(sampled)

    random.shuffle(combined_records)

    output_lines = [line for record in combined_records for line in record]
    output_path = os.path.join(output_dir, f"combined_{i}.fastq")

    with open(output_path, 'w') as out_f:
        out_f.writelines(output_lines)

    print(f"Generated {output_path} with {len(combined_records)} sequences")
