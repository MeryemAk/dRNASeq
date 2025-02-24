
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
  1. A fastq directory containing one subdirectory, which holds the FASTQ files (either in .fastq.gz or .fastq format):
      <img src="images/Tree1.JPG" alt="Beschrijving" width="300" style="display:block; margin-left:auto; margin-right:auto;"/>
  2. A fastq containing multiple subdirectories, each named by a barcode and holding FASTQ files (either in .fastq.gz or .fastq format):
      <img src="images/Tree2.JPG" alt="Beschrijving" width="300" style="display:block; margin-left:auto; margin-right:auto;"/>

## Additional parameters
 * [--barcodes]: (default: none) Comma separated list of barcode numbers that the user wants to analyse if barcodes are present. All barcodes are automatically analysed if barcodes are present but no --barcodes option is provided. Numbers should include the leading 0s. E.g. 03,08,11
 * Parameters related to Quality Control (QC):
    - [--qc]: (default: true) If provided, will perform QC analysis (NanComp and NanoPlot)
    - [--t_qc]: (default: 4) Number of threads used for QC
 * Parameters related to assembly:
   - [--nano_hq] (default: true) Mode for ONT Guppy5+ (SUP mode) and Q20 reads (3-5% error rate)
   - [--gsize]: (default: none) Expected genome size
   - [--meta]: (default: false) Metagenome / Uneven coverage
   - [--asm_coverage]: (default: false) reduced coverage for initial disjointig assembly
   - [--t_assembly]: (default: 4) Number of threads per barcode 
* Parameters related to mapping:
   - [--mapping]: (default: true) If provided, will map sort and index original reads against the assembly
   - [--t_mapping]: (default: 4), minimap2 uses t+1 during mapping) Number of threads used for mapping
* Parameters related to polishing:
   - [--polishing]: (default: true) If provided will polish sequences (requires mapping)
   - [--t_polishing]: (default: 4)  Number of threads used for polishing
   - [--model]:(default: auto selection): Model used for Medaka polishing: {pore}_{device}_{caller variant}_{caller version}, normally automatically selected by Medaka
 * Parameters related to annotation:
   - [--annotation]: (default:true) If true, Prokka will be executed
   - [--t_annotation]: (default: 4) Number of threads used for annotation
 * Parameters related to assembly qc:
   - [--assembly_qc]: (default:true) If provided, will calculate BUSCO scores for the assembly
   - [--t_assembly_qc]: (default: 4) Number of threads used for BUSCO

## Config file
In the config file (nextflow.config), the default parameters can be adjusted, f. ex. the threads used for all the different processes.
Next tot the parameter settings, computing resources can also be modified:
```
executor {
    name = 'local'
    queueSize = 5
    memory = '64 GB'
    cpus = 32
}
```
The maximum memory is limited by default on 64 GB, but should be adjusted to personal computer characteristics.
queueSize is limited to 5, which means that only 5 processes can be executed at once. If the default number of threads are used, this means that 20 (5x4) threads at once will be used. 
The user can also replace queueSize with "cpus" if wanted. 

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


Along with these outputdirectories, 3 Nextflow reports are also generated:
- report: metrics about the workflow execution
- trace: information about each process
- timeline: timeline for all processes

## Other remarks
Each time Nextflow is executed, directories within the work directory are created where the processes run. Don't forget to empty this work direcotry regulary.

