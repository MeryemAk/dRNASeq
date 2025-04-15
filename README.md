# dRNASeq Pipeline  
This repository contains scripts to perform dual RNA-seq analysis. The pipeline is designed for Linux systems and uses Conda environments for reproducibility.

## 📁 Repository Structure  
Scripts are numbered in the order of execution:  
1. `1_linuxsetup.sh` – Sets up Miniconda and creates the required Conda environment.  
    Can be executed manually by following the `manual_setup_linux.txt` document.
2. `2.1_SRR_download.py` - Can be used to automatically download samples from the Sequence Read Archive (SRA) (optional)
3. `2_fastqc.py` – Runs quality control checks
4. ...

## 🛠️ Installation  
Run the setup script to install Miniconda and initialize the environment:

```bash
./linux_setup.sh -h   # View help message first
./linux_setup.sh      # Run setup
```

### 📋 Download samples from SRA (optional)  
The script allows the user to automatically download all samples in fastq.gz format with the fastq-dump tool from the sra-toolkit package.  
The script uses 3 command line arguments:  
1. Path to the file containing SRR accession numbers (one per line).
2. Specify 'SE' for single-end reads or 'PE' for paired-end reads.
3. . Path to the output directory

```bash
python3 2.1_SRR_download.py <accession_file> <SE|PE> [work_dir]
```

## 🧪 FastQC  
Run the `2_fastqc.py` script to perform quality control on all sequences  

```bash
python 2_fastqc.py
```

Notes:
- Place your `.fastq.gz` files (aka sequences) in a folder named `sequences` (or update the folder name in the script).
- The output will be saved in a folder named `QualityControl`.
