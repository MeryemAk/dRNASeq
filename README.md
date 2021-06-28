
# bacterial_assembly
Nextflow script for assembly of bacterial genomes from Nanopore data

## Tools used
* Guppy basecaller (nanopore community)
* Nanocomp https://github.com/wdecoster/nanocomp
* Nanoplot https://github.com/wdecoster/NanoPlot
* Flye https://github.com/fenderglass/Flye
* Minimap2 https://github.com/lh3/minimap2
* Samtools http://www.htslib.org/
* Racon https://github.com/isovic/racon
* Medaka https://github.com/nanoporetech/medaka
* Prokka https://github.com/tseemann/prokka

## General workflow
There are 2 directories that needed to be specified:
- in_dir: the input directory that contains the data that needed to be analysed
- out_dir: the output directory that will contain the results

3 types of input are possible:
1. FAST5 files
2. FASTQ files (multiple FASTQ files that are not merged yet)
3. Merged FASTQ file (one FASTQ file per barcode/sample)

In_dir directory structure for the 3 types of input:
2. FASTQ files: The in_dir must contain a folder named "basecalled" with the FASTQ files. If barcodes are used, a folder for each barcode that contains all the FASTQ-file is expected.
![alt text](https://user-images.githubusercontent.com/56390957/123658980-823e0880-d832-11eb-93bd-eb637d10c8a2.png)

4. Merged FASTQ file: The in_dir must contain a folder named "basecalled" with the merged FASTQ file(s)





