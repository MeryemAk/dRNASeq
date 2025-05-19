## Steps to create the dRNAseq Conda environment  
1. Create Conda environment
``` bash
conda create -n dRNAseq
conda activate dRNAseq
conda config --show channels       # Should show bioconda, conda-forge and defaults
```
2. Install tools
``` bash
conda install -c bioconda nanoplot                                                        # Quality control
pip install NanoComp                                                                      # Quality control
conda install -c nanoporetech -c conda-forge -c bioconda "nanoporetech::pychopper"        # Trimming
conda install bioconda::minimap2                                                          # Mapping
conda install -c bioconda samtools                                                        # Sorting and indexing
conda install -c bioconda subread                                                         # Counting
```  
! error with samtools? --> downgrade python to version 3.10 with following command:
```bash
conda install python=3.10
conda install -c bioconda samtools    # Redownload Samtools
```  
3. Check if all tools installed through
```bash
conda list | grep toolname
```
