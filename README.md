
# Bacterial_assembly
Nextflow script for assembly, polishing, mapping, QC and annotation of bacterial genomes from Nanopore data. The workflow can be used on Linux or Windows. All the different tools run in different containers inside a Nextflow container (optional), which makes the bacterial assembly workflow independent of the platform used. The container images used in the workflow are freely avalaible on the StaPH-B repository (https://hub.docker.com/u/staphb) or are custom made (https://hub.docker.com/u/bikc).

## Tools used in the workflow

![image](https://user-images.githubusercontent.com/56390957/132333064-028d4052-c7d8-4c83-a90c-e12a76356a6e.png)

* Nextflow https://www.nextflow.io/
* Nanocomp https://github.com/wdecoster/nanocomp
* Nanoplot https://github.com/wdecoster/NanoPlot
* Flye https://github.com/fenderglass/Flye
* Minimap2 https://github.com/lh3/minimap2
* Samtools http://www.htslib.org/
* Medaka https://github.com/nanoporetech/medaka
* BUSCO https://busco.ezlab.org/
* Prokka https://github.com/tseemann/prokka

## Possibilities
- QC control of the reads
- Assembly creation
- Mapping of original reads against the assembly
- Polishing of assembly
- Annotation of the assembly
- QC control of the assembly

## Installation Linux
### Prerequisites
On Linux, only Docker is needed: the workflow is started from a Nextflow container. Users can also opt to install Nextflow (https://www.nextflow.io/docs/latest/getstarted.html).

If you want to start the pipeline from FAST5 files, an account on Oxford Nanopore Technologies (ONT) is necessary to use Guppy. Guppy can then be easily implemented in the pipeline if the Docker image of Guppy is provided in de config file. Due to the terms and conditions of ONT, we are not allowed to redistribute Oxford Nanopore software. 

## Quick start
1) Ensure that the assembly.nf script and the nextflow.config file are stored in a folder within your home directory. The data to be analyzed must also be placed in this folder. For details on the required input folder structure, refer to the mandatory parameters. If you don’t necessarily want to use Docker-in-Docker but still want to run Nextflow outside of a Docker container, you can skip this step and proceed to step 3
2) Before running the workflow for the first time, pull the Nextflow Docker image from Docker Hub. If you have already downloaded the image, you can skip this step and proceed to step 3.
```
$ docker pull nextflow/nextflow:21.04.3
```
3) Start the nextflow container. Don't forget to replace "nameoffolder" with your own folder name.
```
$ docker run -it --workdir $PWD -v /var/run/docker.sock:/var/run/docker.sock -v $HOME/"nameoffolder":$HOME/"nameoffolder" nextflow/nextflow:21.04.3 /bin/bash 
```
3) You are now inside the Nextflow container, and the workflow is ready to be executed. Each tool runs in its own separate container, which will be automatically pulled if it is not already available locally when the corresponding process in Nextflow starts.

## Usage
```
nextflow run assembly.nf --in_dir PATH --out_dir PATH
                         [--barcodes]
                         [--qc][--t_qc]
                         [--nano_hq][--gsize][--meta]][--asm_coverage][--t_assembly]
                         [--mapping][--t_mapping]
                         [--polishing][--model][--t_polishing]
                         [--fast_annotation][--t_annotation]
                         [--assembly_qc][--lineage][--busco_path][t_assembly_qc]
                         [--help]
 
For help: nextflow run assembly.nf --help
```

### Mandatory parameters
Two parameters are mandatory:
- in_dir: the input directory that contains the data that needed to be analysed
- out_dir: the output directory that will contain the results

Input should look like one of these examples:
  1. FAST5 files
  2. FASTQ files (multiple FASTQ files that are not merged yet)
  3. A merged FASTQ file (one big FASTQ file per barcode/sample)

For the 3 types of input, the structure must be the following: <br>
  1. FAST5 files: The in_dir must contain a folder named "fast5" with the FAST5 files. 
  2. FASTQ files: The in_dir must contain a folder named "fastq" with the FASTQ files. If barcodes are used, a folder for each barcode that contains all the FASTQ-files for that barcode is expected.
  
  Example of a possible in_dir "bacterial_assembly" with the "fastq" directory that contains FASTQ files for barcode 14 and barcode 16:<br>
  <img src="https://user-images.githubusercontent.com/56390957/132334926-20a5a757-343b-427c-81ef-40f69505a57e.png">

  3. Merged FASTQ file: The in_dir must contain a folder named "merged" inside the folder "fastq" with the merged FASTQ file(s) per barcode(s)
  
  Example of a possible in_dir "bacterial_assembly" with the "fastq" and "merged" directories with the merged FASTQ files of barcode 14 and barcode 16:<br>
  <img src="https://user-images.githubusercontent.com/56390957/132337894-907ac818-f48f-4b2c-8f7d-16faf1879c30.png">

## Additional parameters
 * [--merged]: (default: false) If provided, will expect a fastq/merged direcotry with the merged FASTQ file(s) per barcode(s)
 * [--barcodes]: (default: none) Comma separated list of barcode numbers that the user wants to analyse if barcodes are present. All barcodes are automatically analysed if barcodes are present but no --barcodes option is provided. Numbers should include the leading 0s. E.g. 03,08,11
 * [--basecalling]: (default: false) If provided, this will basecall fast5 files from in_dir
   - [--barcode_kits]: (default: false) If provided, this will demultiplex the samples using the provided kit (f.ex: EXP-NBD196")
   - [--bc_config]: (default: dna_r9.4.1_450bps_fast.cfg) If provided, will use a different config file then default
   - [--skip_qc]: (default: false) If provided, this will not split basecalled reads into pass or fail
   - [--num_callers] (default: 1): Number of callers to use for basecalling (1=4 threads!)
 * [--nanocomp]: (default: true) If provided, will perform nanocomp analysis
 * [--nanoplot]: (default: true) If provided, will perform nanoplot analysis
 * [--t_qc]: (default: 4) Threads used for NanoComp and/or Nanoplot analysis
 * [--flye_assembly]: (default:true) If provided, this will assemble the genomes using Flye
   - [--t_assembly]: (default: 4) Number of threads per barcode 
 * If the flye_assemble option is provided: some extra Flye parameters that can be submitted:
   - [--gsize]: (default: false) Expected genome size
   - [--meta]: (default: false) Metagenome / Uneven coverage
   - [--plasmids]: (default: false) rescue short unassembled plasmids
   - [--asm_coverage]: (default: false) reduced coverage for initial disjointig assembly
 * [--miniasm_assembly]: (default:false) If provided, this will assemble the genomes using Miniasm
 * [--mapping]: (default: true) If provided, will map sort and index original reads against the assembly
   - [--t_mapping]: (default:4), minimap2 uses t+1 during mapping) Number of threads used for mapping
 * [--polishing]: (default:4) If provided will polish sequences (requires mapping)
   - [--t_polishing]: Number of threads used for polishing
   - [--model]:(default: r941_min_fast_g303): Model used for Medaka polishing: {pore}_{device}_{caller variant}_{caller version}
 * [--annotation]: (default:true) If provided, will anotate sequences
   - [--t_annotation]: (default: 4) Number of threads used for annotation

## Config file
In the config file (nextflow.config), the default parameters can be adjusted, f. ex. the threads used for all the different processes.
Next tot the parameter settings, computing resources can also be modified:
```
executor {
    queueSize = 5
    memory = '64 GB'
}
```
The maximum memory is limited by default on 64 GB, but should be adjusted to personal computer characteristics.
queueSize is limited to 5, which means that only 5 processes can be executed at once. If the default number of threads are used, this means that 20 (5x4) threads at once will be used. 
The user can also replace queueSize with "cpus" if wanted. But keep in mind that this constraint concerns only the amount of logical cpus! Almost all the processes in the workflow are multithreaded. If cpus = 24 is used with the default settings, 96(24x4) threads can be used at once!

In the config file the containers can be easily adjusted if other containers and/or version of a tool want to be used.

The report, trace and timeline section in the config file generates reports add the end of the workflow.

## Output
The output is structured in the following way:
1. basecalled (merged FASTQ files)
2. QC (NanoPlot and/or NanoComp output)
3. Assembly (Flye or Miniasm output)
4. Mapping (files necessary for visualisation in IGV)
5. Polishing (
6. Annotation (Output from Prokka)
7. Report (small report with an overview of the tools used and some basic statics (only generated when NanoPlot is used).

Along with these outputdirectories, 3 Nextflow reports are also generated:
- report: metrics about the workflow execution
- trace: information about each process
- timeline: timeline for all processes

## Other remarks
Each time Nextflow is executed, directories within the work directory are created where the processes run. Don't forget to empty this work direcotry regulary.

