# dRNASeq Pipeline  
This repository contains scripts to perform dual RNA-seq analysis. The pipeline is designed for Linux systems and uses Conda environments for reproducibility.

## 📁 Repository Structure  
Scripts are numbered in the order of execution:  
1. `1_linuxsetup.sh` – Sets up Miniconda and creates the required Conda environment.  
    Can be executed manually by following the `manual_setup_linux.txt` document.
2. `2_qc.py` – Runs quality control checks, generates summarizing report and cleans up unnecessary files
3. `3_trimming.py` - Trims long- and short-reads  
4. The mapping consists out of multiple steps  
    4.1. `4.1_create_index.sh` - Downloads bacterial reference genomes and creates an index using `Minimap2`

## 🛠️ 1. Installation  
Run the setup script to install Miniconda and initialize the environment:

```bash
./1_linuxsetup.sh -h   # View help message first
./1_linuxsetup.sh      # Run setup
```  

## 🧪 2. Quality Control (FastQC and LongQC)  
Run the `2_qc.py` script to perform quality control on all sequences. The script will:
1. Filter long reads using either `LongQC` or `FastQC` (change within script!).
    - LongQC for long-reds (ONT)
    - FastQC for short-reads (illumina)
2. Generate a MultiQC report summarizing the quality metrics for all samples.
3. Clean up unnecessary `.zip` files to save disk space.


```bash
python 2_qc.py
```

Notes:
- Input Files: Place your .fastq.gz files (aka sequences) in a folder named sequences.
- Output Folder: The results will be saved in a folder named `QualityControl`. This folder will contain:
    - Individual FastQC reports for each sequence file.
    - A consolidated MultiQC report summarizing the quality metrics for all samples.  

Expected Output Structure:
```bash
QualityControl/
├── sample1.fastqc.html
├── sample2.fastqc.html
├── ...
├── multiqc_report.html
```  

## 3. Trimming  
Run the 3_trimming.py script to trim long-read sequencing data using `Filtlong`. The script applies filters such as minimum read length, percentage of best reads to keep, and a target number of bases. The trimmed reads are compressed into .gz files and saved in the output folder called `trimming`.

```bash
python 3_trimming.py
```

Expected output structure:
```bash
trimming/
├── sample1.fastq.gz
├── sample2.fastq.gz
├── ...
```

## 4 Mapping  
### 4.1 Create index  
Run the 4.1_create_index.sh script to download bacterial genomes from the NCBI RefSeq database and create an index using `Minimap2`.

```bash
./4.1_create_index.sh
```

Notes:
- Add or delete acceccion numbers within the script depending on which bacterial reference genomes you want to download and create an index of.
- All reference genomes will be concatenated in `all_genomes.fna` which will then be used to create the `bacterial_index.mmi` index with `Minimap2`