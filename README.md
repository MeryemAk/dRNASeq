# dRNASeq Pipeline  
This repository contains scripts to perform dual RNA-seq analysis. The pipeline is designed for Linux systems and uses Conda environments for reproducibility.

## 📁 Repository Structure  
Scripts are numbered in the order of execution:  
1. `1_linuxsetup.sh` – Sets up Miniconda and creates the required Conda environment.  
    Can be executed manually by following the `manual_setup_linux.txt` document.
2. `2.1_SRR_download.py` - Can be used to automatically download samples from the Sequence Read Archive (SRA) (optional)
3. `2_qc.py` – Runs quality control checks, generates summarizing report and cleans up unnecessary files
4. ...

## 🛠️ Installation  
Run the setup script to install Miniconda and initialize the environment:

```bash
./1_linuxsetup.sh -h   # View help message first
./1_linuxsetup.sh      # Run setup
```

### 📋 Download samples from SRA (optional)  
The script allows the user to automatically download all samples in `fastq.gz` format using the `fastq-dump` tool from the `sra-toolkit` package.  
The script uses 3 command line arguments:  
1. Path to the file containing SRR accession numbers (one per line).
2. Specify 'SE' for single-end reads or 'PE' for paired-end reads.
3. Path to the output directory  

```bash
python3 2.1_SRR_download.py <accession_file> <SE|PE> [work_dir]
```

## 🧪 Quality Control (FastQC and Filtlong)  
Run the `2_qc.py` script to perform quality control on all sequences. The script will:
1. Filter long reads using either `Filtlong` or `FastQC` (change within script!).
2. Generate a MultiQC report summarizing the quality metrics for all samples.
3. Clean up unnecessary `.zip` files to save disk space.


```bash
python 2_qc.py
```

Notes:
- Input Files: Place your .fastq.gz files (aka sequences) in a folder named sequences.
- Output Folder: The results will be saved in a folder named QualityControl. This folder will contain:
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
