
# Bacterial_assembly
Nextflow script for assembly of bacterial genomes from Nanopore data. This workflow can be used on Linux or Windows. All the different tools run in different containers, which makes the bacterial assembly workflow independent of the platform used. 

## Tools used in the workflow

![image](https://user-images.githubusercontent.com/56390957/132333064-028d4052-c7d8-4c83-a90c-e12a76356a6e.png)

* Nextflow https://www.nextflow.io/
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

## Installation
### Prerequisites
On Windows, only Docker is needed for the installation of this tool (https://docs.docker.com/get-docker/).
On Linux, Docker is needed. One may choose to install Nextflow (https://www.nextflow.io/docs/latest/getstarted.html) or use the Nextflow container.

## Quick start
1) Make sure that the assembly.nf script and the nextflow.config file are in your home directory. The data that will be analysed also needs to be in this folder (for more details on de structure of the input folder, see mandatory parameters).
2) Pull the Docker image of nextflow from Dockerhub (This is only necesssary the first time the workflow is used):
```
$ docker pull nextflow:nextflow
```
3) Start the nextflow container:
```
$ docker run -it --workdir $PWD -v /var/run/docker.sock:/var/run/docker.sock -v $HOME:$HOME  nextflow/nextflow /bin/bash 
```
3) Now you are in the Nextflow container and the workflow can be executed. All the different tools run in different containers. These will be pulled automatically (when you don't have them yet locally) when the corresponding process in Nextflow is started.

## Usage
```
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
```

### Mandatory parameters
Two parameters are mandatory:
- in_dir: the input directory that contains the data that needed to be analysed
- out_dir: the output directory that will contain the results

3 types of input are possible:
  1. FAST5 files
  2. FASTQ files (multiple FASTQ files that are not merged yet)
  3. A merged FASTQ file (one big FASTQ file per barcode/sample)

For the 3 types of input, the structure must be the following: <br>
  1. <br>
  2. FASTQ files: The in_dir must contain a folder named "basecalled" with the FASTQ files. If barcodes are used, a folder for each barcode that contains all the FASTQ-files for that barcode is expected.
  
  Example of a possible in_dir "bacterial_assembly" with the "basecalled" directory that contains FASTQ files for barcode 14 and barcode 16:<br>
  <img src="https://user-images.githubusercontent.com/56390957/132334926-20a5a757-343b-427c-81ef-40f69505a57e.png">

  3. Merged FASTQ file: The in_dir must contain a folder named "merged" with the folder "basecalled" with the merged FASTQ file(s) per barcode(s)
  
  Example of a possible in_dir "bacterial_assembly" with the "basecalled" and "merged" directory with the merged FASTQ files of barcode 14 and barcode 16:<br>
  <img src="https://user-images.githubusercontent.com/56390957/132337894-907ac818-f48f-4b2c-8f7d-16faf1879c30.png">

## Additional parameters
 * [--merged]: (default:false) If provided, will expect a basecalled/merged direcotry with the merged FASTQ file(s) per barcode(s)
 * [--barcodes]: (default:none) Comma separated list of barcode numbers that are expected, if barcodes are provided. Numbers should include the leading 0s. E.g. 03,08,11. 
 * [--nanocomp]: (default:true) If provided, will perform nanocomp analysis
 * [--nanoplot]: (default:true) If provided, will perform nanoplot analysis
 * [--assemble]: (default:true) If provided, this will assemble the genomes using Flye
   - [--t_assemble]: (default: 8) Number of threads per barcode 
 * If the assemble option is provided: some extra Flye parameters that can be submitted:
   - [--gsize]: Expected genome size
   - [--meta]: Metagenome / Uneven coverage
   - [--plasmids]: rescue short unassembled plasmids
   - [--asm_coverage]: reduced coverage for initial disjointig assembly
 * [--mapping]: (default: true) If provided, will map sort and index original reads against the assembly
   - [--t_mapping]: (default: 3, minimap2 uses t+1 during mapping) Number of threads used for mapping
 * [--polishing]: (default:4) If provided will polish sequences (requires mapping)
   - [--t_polishing]: Number of threads used for polishing
   - [--model]:(default: r941_min_fast_g303): Model used for Medaka polishing: {pore}_{device}_{caller variant}_{caller version}
 * [--annotation]: (default:true) If provided, will anotate sequences
   - [--t_polishing]: (default: 4) Number of threads used for annotation


