# Dual RNA-Seq workflow
This Nextflow pipeline is designed for the analysis of dual RNA-seq data from vaginal swabs. It performs preprocessing, quality control (QC), mapping quantification and taxonomic classification of reads from host (Homo sapiens), yeast (Candida albicans) and bacterial genomes.

## Tools used in the workflow

<img src="images/Pipeline.png" alt="Pipeline" width="800" style="display:block; margin-left:auto; margin-right:auto;"/>

* [Nextflow](https://www.nextflow.io/)
* [Nanopack](https://github.com/wdecoster/nanopack)
* [Pychopper](https://github.com/epi2me-labs/pychopper)
* [Minimap2](https://github.com/lh3/minimap2)
* [MMseqs2](https://github.com/soedinglab/MMseqs2)
* [gffread](https://github.com/gpertea/gffread)
* [Samtools](http://www.htslib.org/)
* [Kraken2](https://ccb.jhu.edu/software/kraken/)
* [bambu](https://hub.docker.com/r/mathiasverbeke/bambu_runner)

## Possibilities
- Filter 16S/18S and 23S/28S rRNA from samples
- Trim reads
- QC control of the reads
- Mapping of reads against Homo sapiens, Candida albicans and bacterial reference genomes
- Quantification
- Taxonomic classification of unmapped reads

## Wiki pages
The workflow can be executed in two ways:  
- Using the singular scripts
- Using Nextflow (not available yet)

For more information on how to start, refer to the Standard Operating Procedures of [Nextflow](https://github.com/MeryemAk/dRNASeq/wiki/Standard-Operating-Procedure-for-Nextflow) or the [scripts](https://github.com/MeryemAk/dRNASeq/wiki/Standard-Operating-Procedure-for-scripts).

## Repository structure
* `1.data/`: store data in this directory in fastq format  
* `images/`: images used in documentation  
* `misc/`: miscellaneous files (other documenatation relevant to the GitHub repository)  
* `nanocomp/`: Docker file to build nanocomp  
* `nanopack/`: Docker file to build nanopack  
* `reference_genomes/`: store reference genomes and annotation files here. Also contains the `create_bacteria_index.sh` script to download bacterial files  
* `scripts/`: contains all script needed to run the pipeline individually  
* `environment.yml`: file used to build the Conda environment  
* `main.nf`: main nextflow script that runs all steps  
* `nextflow.config`: configuration file for main.nf  
* `README.md`: this document  

## Other remarks
Each time Nextflow is executed, directories within the work directory are created where the processes run. Don't forget to empty this work direcotry regulary.

## Questions? meryem191101@gmail.com
