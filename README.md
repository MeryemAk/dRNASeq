
# bacterial_assembly
Nextflow script for assembly of bacterial genomes from nanopore data

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
3 types of input are possible.
- FAST5 files
- FASTQ files (mutliple FASTQ files per barcode that are not merged yet)
- Merged FASTQ file (one FASTQ file per barcode)



