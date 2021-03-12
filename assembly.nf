#!/usr/bin/env nextflow 

/* Help message*/
if(params.help) {
    print helpMessage()
    exit 0
}
else {
    print parameterShow()
}

/* parameters validation */
error=false

if (!params.in_dir){
    println "${c_red}params.in_dir was empty - no input directory specified${c_reset}"
    error=true
} else if (!file(params.in_dir).isDirectory()){
    println "${c_red}params.in_dir is not a directory - no input directory specified${c_reset}"
    error=true
}

if (!params.out_dir){
    println "${c_red}params.out_dir was empty - no output directory specified${c_reset}"
    error=true
} else if (!file(params.out_dir).isDirectory()){
    println "${c_red}params.out_dir is not a directory - no output directory specified${c_reset}"
    error=true
}

if (error) {
    exit 0
}

/* processes */
process guppy {

    container 'genomicpariscentre/guppy:4.4.2'
    
    publishDir "${params.out_dir}/basecalled", mode: "copy"
    
    input:
    path x from Channel.fromPath(params.in_dir+"/fast5_all")
    
    output:
    path "merged/*.fastq" into fastq_w_guppy
    path "*" into basecall_result

    when:
    params.basecall

    script:
    barcode_kits = params.barcode_kits ? "--barcode_kits $params.barcode_kits" : ""
    skip_qc = params.skip_qc ? "" : "--qscore_filtering"

    """
    guppy_basecaller \
    --input_path ${params.in_dir}/fast5_all \
    --save_path ./ \
    -c ${params.bc_config} \
    --num_callers ${params.num_callers} \
    $barcode_kits \
    $skip_qc 

    # when qc is skipped --> no pass directory
    mkdir merged
    if [[ ${params.skip_qc} == false ]]; then
        dir="./pass/*"
    else
        dir="./*"
    fi

    for d in \$dir
    do
        # merge without barcodes
        if [[ -f \$d && \$d == *.fastq ]]; then
            cat \$dir*.fastq > ./merged/sample.fastq
            break
        # merge with barcodes
        elif [[ \$d =~ barcode* ]]; then
            cat \$d"/"* > ./merged/\$(basename \$d)".fastq"
        fi
    done
    """
}

process merge {

    publishDir "${params.out_dir}/basecalled/", mode: "copy"
    
    input: 
    path x from Channel.fromPath(params.in_dir+"/basecalled").collect()

    when:
    !params.no_merge && !params.basecall

    output:
    path "merged/*.fastq" into fastq_wo_guppy

    """
    mkdir merged
    for d in ./basecalled/*
    do
        if [[ -f \$d ]]; then
            cat ./basecalled/*.fastq > ./merged/sample.fastq
            break
        elif [[ \$d =~ barcode* ]]; then
            cat \$d"/"* > ./merged/\$(basename \$d)".fastq"
        fi
    done
    """
}

if (params.no_merge && !params.basecall) {
	fastq_merged = Channel.fromPath(params.in_dir+"/basecalled/merged/*.fastq")
	fastq_w_guppy = Channel.empty()
	fastq_wo_guppy = Channel.empty()
}
else {
	fastq_merged = Channel.empty()
}

fastqs = fastq_w_guppy.mix(fastq_wo_guppy, fastq_merged)
fastqs.into{fastq_nanoplot;fastq_nanocomp;fastq_assembly}

process nanoplot {

    container 'staphb/nanoplot:1.33.0'

    publishDir "${params.out_dir}/qc/", mode:"copy"

    input:
    path fq from fastq_nanoplot

    output:
    path "*" into nanoplot_result

    when:
    params.nanoplot

    script:
    """
    NanoPlot -t 4 -o ./nanoplot --fastq_rich $fq
    """
}

process nanocomp {

    container 'nanocomp:latest'

    publishDir "${params.out_dir}/qc/", mode:"copy"

    input:
    path fq from fastq_nanocomp

    output:
    path "*" into nanocomp_result
    
    when:
    params.nanocomp 

    script:
    """
    NanoComp -t 4 -o ./nanocomp -f png --plot violin --fastq $fq
    """

}

process flye_assembly {

    container 'staphb/flye:2.8'

	publishDir "${params.out_dir}/assembly/", mode: "copy"

	input: 
	path fq from fastq_assembly.flatten()

	output:
	path "*" optional true into assembly
    tuple file(fq), file("*/assembly.fasta") optional true into (flye_prokka, mapping)

	when:
	params.assemble

    script:

	gsize = params.gsize ? "--genome-size $params.gsize" : ""

    """
    IFS=','
    read -ra BCS <<<"${params.barcodes}"
    # flye zonder barcodes
    if [[ \${#BCS[@]} -eq 0 ]]; then
        flye --nano-raw $fq $gsize --out-dir ./${fq.simpleName} -t ${params.assemblyP} -i 5
    #flye met barcodes
    else
        count=0
        for i in \${BCS[@]}
        do
            BCS[\$count]="barcode"\$i".fastq"
            echo "barcode"\$i".fastq"
            count=\$((count + 1))
        done
        if [[ " \${BCS[@]} " =~ " $fq " ]]; then
            flye --nano-raw $fq $gsize --out-dir ./${fq.simpleName} -t ${params.assemblyP} -i 5
        fi

    fi
    """
} 

process mapping {

    container 'staphb/minimap2:2.17'

    publishDir "${params.out_dir}/mapping/", mode: "copy"

    input:
    tuple file(fq), file(assembly) from mapping

    output:
    path "*" into sam_files
    tuple file(fq), file(assembly), file("*.sam") into racon

    when:
	params.mapping
    
    """
    minimap2 -a $assembly $fq > ${fq.simpleName}.sam
    """
}

process sam_to_bam {

    container 'staphb/samtools:1.11'

    publishDir "${params.out_dir}/mapping", mode: "copy"

    input:
    path x from sam_files

    output:
    path "*" into bam_files

    when:
	params.mapping
    
    """
    samtools view -S -b $x > ${x.simpleName}.bam
    """
}

process sort_bam {

    container 'staphb/samtools:1.11'

    publishDir "${params.out_dir}/mapping", mode: "copy"

    input:
    path x from bam_files

    output:
    path "*" into sorted_bam_files

    when:
	params.mapping
    
    """
    samtools sort $x -o ${x.simpleName}_sorted.bam
    """
}

process index_bam {
    
    container 'staphb/samtools:1.11'

    publishDir "${params.out_dir}/mapping", mode: "copy"

    input:
    path x from sorted_bam_files

    output:
    path "*" into indexed_bam_files

    when:
	params.mapping
    
    """
    samtools index $x
    """
}

process racon {

    container 'staphb/racon:1.4.20'

    publishDir "${params.out_dir}/polishing/racon", mode: "copy"

    input:
    tuple file(fq), file(assembly), file(alignment) from racon

    output:
    path "*" into racon_result
    tuple file(fq), file("*_racon.fasta") into medaka

    when:
    params.polishing

    """
    # -m 8 -x -6 -g -8 -w 500: recommended parameters before using medaka
    racon $fq $alignment $assembly -t ${params.threads_polishing} -m 8 -x -6 -g -8 -w 500 > ${fq.simpleName}_racon.fasta 
    """

}

process medaka {

    container 'staphb/medaka:1.2.0'

    publishDir "${params.out_dir}/polishing/medaka", mode: "copy"

    input:
    tuple file(fq), file(assembly) from medaka

    output:
    path "*" into medaka_result
    tuple file(fq), file ("*/consensus.fasta") into medaka_prokka

    when:
    params.polishing

    """
    medaka_consensus -i $fq -d $assembly -o ./${fq.simpleName} -t ${params.threads_polishing} -m ${params.model} 
    """

}

if (params.polishing) {
    prokka_tuple = medaka_prokka
}
else {
    prokka_tuple = flye_prokka
}

process prokka_annotation {

    container 'staphb/prokka:1.14.5'
    
    publishDir "${params.out_dir}/prokka/", mode: "copy"

    input:
    tuple file(fq), file(assembly) from prokka_tuple

	output:
	path "*" into prokka_result

	when:
	params.annotation

    """
    prokka --outdir ./${fq.simpleName} --prefix ${fq.simpleName} $assembly
    """
}

/* help */

def helpMessage() {

    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_orange = params.monochrome_logs ? '' : "\033[0;202m";
    c_gray = params.monochrome_logs ? '' : "\033[0;241m";


    return """\
    ==================================================
    Parameters:

    ${c_cyan}General:${c_reset}                        | ${c_blue}Must supply parameters${c_reset}
    • in_dir:                       | ${c_yellow}Input directory (directory with fast5_all and/or fastq_pass folders)${c_reset}
    • out_dir:                      | ${c_yellow}Output directory with results${c_reset}

    ${c_cyan}Basecalling:${c_reset}                    | ${c_blue}Only when providing --basecall option${c_reset}
    • basecall:                     | ${c_yellow}If provided, this will basecall fast5 files from in_dir ${c_reset}
    • barcode_kits:                 | ${c_yellow}If provided, this will demultiplex the samples using the provided kit${c_reset}
    • bc_config:                    | ${c_yellow}If provided, will use a different config file then default${c_reset}
        default: dna_r9.4.1_450bps_fast.cfg
    • skip_qc:                      | ${c_yellow}If provided, this will not split basecalled reads into pass or fail${c_reset}
    • num_callers (default:1):      | ${c_yellow}Number of callers to use for basecalling (1=4 threads!)${c_reset}

    ${c_cyan}Merging:${c_reset}  
    • no_merge:                     | ${c_yellow}If provided: will expect a /basecalled/ folder with fastq files for each barcode / a sample${c_reset}
                                    | ${c_yellow}If not provided: will expect a /basecalled/merged/ folder with merged fastq files for each barcode${c_reset}
                                    | ${c_yellow}/ a sample${c_reset}
   
    ${c_cyan}QC:${c_reset}                             | ${c_blue}Quality Control steps${c_reset}
    • nanocomp:                     | ${c_yellow}If provided, will perform nanocomp analysis ${c_reset}
    • nanoplot:                     | ${c_yellow}If provided, will perform nanoplot analysis ${c_reset}

    ${c_cyan}Assembly:${c_reset}                       | ${c_blue}Only when providing --assemble option${c_reset}
    • assemble:                     | ${c_yellow}If provided, this will assemble the genomes using Flye ${c_reset}
    • gsize:                        | ${c_yellow}Expected genome size (not mandatory) ${c_reset}
    • assemblyP (default: 8):       | ${c_yellow}Number of threads per barcode (use max: 32/nbarcodes) ${c_reset}
    • barcodes:                     | ${c_yellow}Comma separated list of barcode numbers that are expected, if barcodes are provided${c_reset}
                                    | ${c_yellow}Numbers should include the leading 0s. E.g. 03,08,11${c_reset}

    ${c_cyan}mapping:${c_reset}                        | ${c_blue}If provided, will map, sort and index sequences ${c_reset}

    ${c_cyan}polishing:${c_reset}                      | ${c_blue}If provided, will polish sequences (requires mapping) ${c_reset}
    • threads_polishing             | ${c_yellow}Number of threads used for polishing ${c_reset}
    • model:                        | ${c_yellow}Model used for Medaka polishing: {pore}_{device}_{caller variant}_{caller version} ${c_reset}
        default: r941_min_fast_g303

    ${c_cyan}annotation:${c_reset}                     | ${c_blue}If provided, will anotate sequences (requires mapping) ${c_reset}
   
    • help                          | ${c_yellow}Show this${c_reset}
                        
    
    Run using: nextflow run assembly.nf --in_dir {Input directory} --out_dir {Output directory}
    ==================================================
    """

}

def parameterShow() {

    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_orange = params.monochrome_logs ? '' : "\033[0;202m";
    c_gray = params.monochrome_logs ? '' : "\033[0;241m";

    text = """\
    ==================================================
    Provided parameters:

    ${c_cyan}General:${c_reset} 
    • in_dir:                       | ${c_yellow}${params.in_dir}${c_reset}
    • out_dir:                      | ${c_yellow}${params.out_dir}${c_reset}

    ${c_cyan}Basecalling:${c_reset} 
    • basecall:                     |"""

    if (params.basecall) {
        text = text + """ ${c_yellow}${params.basecall}${c_reset}
    • barcode_kits:                 | ${c_yellow}${params.barcode_kits}${c_reset}
    • bc_config:                    | ${c_yellow}${params.bc_config}${c_reset}
    • skip_qc:                      | ${c_yellow}${params.skip_qc}${c_reset}
    • num_callers:                  | ${c_yellow}${params.num_callers}${c_reset}

    """
    }
    else {
        text = text+""" ${c_yellow}${params.basecall}${c_reset}

    """
    }

    text = text + """${c_cyan}Merge:${c_reset}
    • no_merge:                     | ${c_yellow}${params.no_merge}${c_reset}

    """

    text = text + """${c_cyan}QC:${c_reset}
    • nanocomp:                     | ${c_yellow}${params.nanocomp}${c_reset}
    • nanoplot:                     | ${c_yellow}${params.nanoplot}${c_reset}

    ${c_cyan}Assembly:${c_reset}
    • assemble:                     |"""

    if (params.assemble) {
        text = text + """ ${c_yellow}${params.assemble}${c_reset}
    • gsize:                        | ${c_yellow}${params.gsize}${c_reset}
    • assemblyP:                    | ${c_yellow}${params.assemblyP} ${c_reset}
    • barcodes:                     | ${c_yellow}${params.barcodes} ${c_reset}

    """
    }
    else {
        text = text+""" ${c_yellow}${params.assemble}${c_reset}

    """
    }

    text = text + """${c_cyan}Mapping:${c_reset}
    • mapping:                      | ${c_yellow}${params.mapping}${c_reset}

    """

    text = text + """${c_cyan}Polishing:${c_reset}
    • polishing:                     |"""

    if (params.polishing) {
        text = text + """ ${c_yellow}${params.polishing}${c_reset}
    • model:                        | ${c_yellow}${params.model} ${c_reset}
    • threads_polishing:            | ${c_yellow}${params.threads_polishing} ${c_reset}

    """
    }
    else {
        text=text+""" ${c_yellow}${params.polishing}${c_reset}

    """
    }

    text = text + """${c_cyan}Annotation:${c_reset}
    • annotation:                   | ${c_yellow}${params.annotation}${c_reset}

    """

    text=text+"""Output directory:  | ${c_yellow}${params.out_dir}${c_reset}
    =================================================="""


    return text;
}

