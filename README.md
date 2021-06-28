
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

## Possibilities
- QC control of the reads
- Assembly creation
- Mapping of original reads against the assembly
- Polishing of assembly
- Annotation of the assembly

## Usage
nextflow run assembly.nf --in_dir PATH --out_dir PATH
                         [--basecall][--barcode_kits][--bc_config][--skip_qc][--num_callers]
                         [--no_merge]
                         [--nanocomp][--nanoplot]
                         [--assemble][--gsize][--meta][--plasmids][--asm_coverage][--assemblyP]
                         [--mapping]
                         [--polishing][--model][--trheads_polishing]
                         [--annotation][--threads_annotation]
                         [--help]
 optional arguments:
  --help 


## Input parameters
There are 2 directories that needed to be specified:
- in_dir: the input directory that contains the data that needed to be analysed
- out_dir: the output directory that will contain the results

3 types of input are possible:
  1. FAST5 files
  2. FASTQ files (multiple FASTQ files that are not merged yet)
  3. Merged FASTQ file (one FASTQ file per barcode/sample)

in_dir directory structure for the 3 types of input: <br>
  1. <br>
  2. FASTQ files: The in_dir must contain a folder named "basecalled" with the FASTQ files. If barcodes are used, a folder for each barcode that contains all the FASTQ-files is expected.
<p align="left" width="100%">
  basecalled directory with the FASTQ files per barcode: <br>
  <img width="20%" src="https://user-images.githubusercontent.com/56390957/123658980-823e0880-d832-11eb-93bd-eb637d10c8a2.png">
  <img width="30%" src="https://user-images.githubusercontent.com/56390957/123661149-9c78e600-d834-11eb-9f3a-0c245b3ce6c8.png">
</p>

  3. Merged FASTQ file: The in_dir must contain a folder named "basecalled" with the merged FASTQ file(s) <br>

## Additional parameters
### 01. Basecalling




 
 





