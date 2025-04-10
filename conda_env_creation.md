## Steps to create the dRNAseq Conda environment  
1. Create Conda environment
``` bash
conda create -n dRNAseq
conda activate dRNAseq
conda config --show channels       # Should show bioconda, conda-forge and defaults
```
2. Install tools
``` bash
conda install -c bioconda fastqc
conda install -c bioconda trimmomatic
conda install -c bioconda #mappingtool?
conda install -c bioconda samtools
conda install -c bioconda #countingtool?
```  
! error with samtools? --> downgrade python to version 3.10 with following command:
```bash
conda install python=3.10
conda install -c bioconda samtools    # Redownload Samtools
```
