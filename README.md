# dRNASeq Pipeline

This repository contains scripts to perform dual RNA-seq analysis. The pipeline is designed for Linux systems and uses Conda environments for reproducibility.

## 📁 Repository Structure

Scripts are numbered in the order of execution:

1. `1_linuxsetup.sh` – Sets up Miniconda and creates the required Conda environment.  
    Can be executed manually by following the `manual_setup_linux.txt` document.
3. `2_fastqc.py` – Runs quality control checks
4. `script_2.py` – ...
5. ...

## 🛠️ Installation

Run the setup script to install Miniconda and initialize the environment:

```bash
./linux_setup.sh -h   # View help message first
./linux_setup.sh      # Run setup

## 🧪 FastQC

Run the `2_fastqc.py` script to perform quality control on all sequences


```bash
python 2_fastqc.py
```

Notes:
- Place your `.fastq.gz` files (aka sequences) in a folder named `sequences` (or update the folder name in the script).
- The output will be saved in a folder named `QualityControl`.