#!/usr/bin/env python3

import os

# PARAMETERS ####
outdir = "./mapping_result"
indir = "/home/meryema/fastq_dump/data"
#"/media/sf_SF/00-stage/documents/test/3.mapping/data"
human_index = "/data/igenomes/Homo_sapiens/GRCh38.p14/minimap2/genome.mmi"
#"/media/sf_SF/00-stage/documents/test/3.mapping/human_refgenome.fna"         # --> GCF_000001405.40
candida_index = "/data/igenomes/Candida_albicans/minimap2/genome.mmi" 
#"/media/sf_SF/00-stage/documents/test/3.mapping/candida_refgenome.fna"     # --> GCF_000182965.3
bacterial_index = "/home/meryema/mapping/bacterial_index.mmi" 
#"/media/sf_SF/00-stage/documents/test/3.mapping/bacterial_index.mmi"     # --> E.coli and S.aureus

# Create output directory if it doesn't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)
    print(f"Created directory: {outdir}")
else:
    print(f"Directory already exists: {outdir}")

# Loop over files from indir
for file in os.listdir(indir):
    sample_name = os.path.splitext(file)[0]  # Get the sample name without extension
    input_file = os.path.join(indir, file)  # Make filename
    
    print(f"===== processing file: {sample_name} =====")

    ################################################################################################
    # Step 1: Map to Human Genome
    print("=== mapping to human reference genome ===")
    human_output = os.path.join(outdir, f"{sample_name}_human.sam")
    if not os.path.exists(human_index):
        print(f"Error: Human refernce genome file not found: {human_index}")
        exit(1)
    else:
        print(f"Human reference genome file found: {human_index}")
        command_human = f"minimap2 -ax splice --secondary=no {human_index} {input_file} > {human_output}"
        os.system(command_human)
        print(f"Mapped {input_file} to human genome: {human_output}")

    print(" ")
    # Extract unmapped and mapped reads from Human Genome
    print("===== extracting unmapped and mapped reads from human genome =====")
    unmapped_human = os.path.join(outdir, f"{sample_name}_human_unmapped.fastq")
    mapped_human = os.path.join(outdir, f"{sample_name}_human_mapped.sam")
    os.system(f"samtools fastq -f 4 {human_output} > {unmapped_human}")  # Unmapped reads
    os.system(f"samtools view -h -F 4 {human_output} > {mapped_human}")      # Mapped reads

    # Sort and index the mapped reads
    sorted_human_bam = os.path.join(outdir, f"{sample_name}_mapped_human_sorted.bam")
    command_sort_human = f"samtools sort {mapped_human} -o {sorted_human_bam}"
    os.system(command_sort_human)
    command_index_human = f"samtools index {sorted_human_bam}"
    os.system(command_index_human)
    print(f"Sorted and indexed human genome mapping: {sorted_human_bam}")

    print(" ")
    ################################################################################################
    # Step 2: Map to Candida Genome
    print("===== mapping to candida reference genome =====")
    # Check if the Candida index file exists
    candida_output = os.path.join(outdir, f"{sample_name}_candida.sam")
    if not os.path.exists(candida_index):
        print(f"Error: Candida reference genome file not found: {candida_index}")
        exit(1)
    else:
        print(f"Candida reference genome file found: {candida_index}")
        command_candida = f"minimap2 -ax splice --secondary=no {candida_index} {unmapped_human} > {candida_output}"
        os.system(command_candida)
        print(f"Mapped unmapped reads to Candida genome: {candida_output}")

    print(" ")
    # Extract unmapped and mapped reads from Candida Genome
    print("===== extracting unmapped and mapped reads from candida genome =====")
    unmapped_candida = os.path.join(outdir, f"{sample_name}_candida_unmapped.fastq")
    mapped_candida = os.path.join(outdir, f"{sample_name}_candida_mapped.sam")
    os.system(f"samtools fastq -f 4 {candida_output} > {unmapped_candida}")  # Unmapped reads
    os.system(f"samtools view -h -F 4 {candida_output} > {mapped_candida}")     # Mapped reads

    print(" ")

    # Sort and index the mapped reads
    sorted_candida_bam = os.path.join(outdir, f"{sample_name}_mapped_candida_sorted.bam")
    command_sort_candida = f"samtools sort {mapped_candida} -o {sorted_candida_bam}"
    os.system(command_sort_candida)
    command_index_candida = f"samtools index {sorted_candida_bam}"
    os.system(command_index_candida)
    print(f"Sorted and indexed Candida genome mapping: {sorted_candida_bam}")

    print(" ")
    ################################################################################################
    # Step 3: Map to Bacterial Index
    print("===== mapping to bacterial reference genome =====")
    bacterial_output = os.path.join(outdir, f"{sample_name}_bacterial.sam")
    if not os.path.exists(bacterial_index):
        print(f"Error: Bacterial index file not found: {bacterial_index}")
        exit(1)
    else:
        print(f"Bacterial index file found: {bacterial_index}")
        command_bacterial = f"minimap2 -ax map-ont --secondary=no {bacterial_index} {unmapped_candida} > {bacterial_output}"
        os.system(command_bacterial)
        print(f"Mapped unmapped reads to bacterial genome: {bacterial_output}")
    
    print(" ")
    # Extract unmapped reads from bacterial index
    # Extract unmapped and mapped reads from Bacterial Genome
    print("===== extracting unmapped and mapped reads from bacterial genome =====")
    unmapped_bacteria = os.path.join(outdir, f"{sample_name}_bacteria_unmapped.fastq")
    mapped_bacteria = os.path.join(outdir, f"{sample_name}_bacteria_mapped.sam")
    os.system(f"samtools fastq -f 4 {bacterial_output} > {unmapped_bacteria}")  # Unmapped reads
    os.system(f"samtools view -h -F 4 {bacterial_output} > {mapped_bacteria}")     # Mapped reads

    print(" ")
    
    # Sort and index the mapped reads
    sorted_bacteria_bam = os.path.join(outdir, f"{sample_name}_mapped_bacterial_sorted.bam")
    command_sort_bacteria = f"samtools sort {mapped_bacteria} -o {sorted_bacteria_bam}"
    os.system(command_sort_bacteria)
    command_index_bacteria = f"samtools index {sorted_bacteria_bam}"
    os.system(command_index_bacteria)
    print(f"Sorted and indexed bacterial genome mapping: {sorted_bacteria_bam}")

    print("="*100)
print("Mapping completed for all files.")