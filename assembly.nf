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
    println "${c_red}\nparams.in_dir was empty - no input directory specified${c_reset}"
    error=true
} 
else if (!file(params.in_dir).isDirectory()){
    println "${c_red}\nparams.in_dir is not a directory - no input directory specified${c_reset}"
    error=true
}

if (!params.out_dir){
    println "${c_red}\nparams.out_dir was empty - no output directory specified${c_reset}"
    error=true
} 
else if (!file(params.out_dir).isDirectory()){
    println "${c_red}\nparams.out_dir is not a directory - no output directory specified${c_reset}"
    error=true
}

if (!params.basecall && !params.no_merge && (!file("${params.in_dir}/basecalled/").isDirectory())) {
    println "${c_red}\nparams.in_dir needs to contain a basecalled directory with fastq files${c_reset}"
    println "${c_red}\nMake a basecalled directory inside the params.in_dir${c_reset}"
    error=true
}


if (!params.basecall && !params.no_merge && (file("${params.in_dir}/basecalled/merged").isDirectory())) {
    println "${c_red}\n params.in_dir contains a basecalled/merged directory, fastqc's are already merged! ${c_reset}"
    println "${c_red}\n Set params.no_merge to true ${c_reset}"
    error=true
}

if (!params.basecall && params.no_merge && (!file("${params.in_dir}/basecalled/merged").isDirectory())) {
    println "${c_red}\n params.in_dir needs to contain a basecalled/merged directory when params.no_merge = true${c_reset}"
    error=true
}
else if (!params.basecall && params.no_merge && (file("${params.in_dir}/basecalled/merged").isDirectory())) {
    count = 0
    file("${params.in_dir}/basecalled/merged").eachFile {
        if (it.name.endsWith('.fastq')) {
            count += 1
        }
    }
    if (count == 0) {
        println "${c_red}\nNo fastq files found in the basecalled/merged directory${c_reset}"
        error=true
    }
}

if (params.asm_coverage && !params.gsize){
    println "${c_red}\nIf params.asm_coverage is specified, params.gsize is mandatory${c_reset}"
    error=true
}
if (error) {
    exit 0
}

/* processes */
process guppy {

    container 'genomicpariscentre/guppy:4.4.2'
    
    publishDir "${params.out_dir}/01.basecalled", mode: "copy"
    
    input:
    path x from Channel.fromPath(params.in_dir+"/fast5_all")
    
    output:
    path "merged/*.fastq" into fastq_w_guppy
    path "*" 

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

    IFS=','
    read -ra BCS <<<"${params.barcodes}"
    # als maar 1 barcode meegegeven wordt met leading 0: octaal voor bash
    if [[ \${#BCS[@]} -eq 1 ]] && [[ \${#BCS[0]} -eq 1 ]]; then
        BCS[0]="0"\${BCS[0]}
    fi
    
    for d in \$dir
    do
        # merge without barcodes
        if [[ -f \$d && \$d == *.fastq ]]; then
            cat \$dir*.fastq > ./merged/sample.fastq
            break
        # merge with barcodes
        elif [[ \$d =~ barcode* ]]; then
            # merge all the barcodes
            if [[ \${#BCS[@]} -eq 0 ]]; then
                cat \$d"/"* > ./merged/\$(basename \$d)".fastq"
            # merge only the barcodes from params.barcodes
            else
                bc=\${d##*barcode}
                if [[ " \${BCS[@]} " =~ " \$bc " ]]; then
                    cat \$d"/"* > ./merged/\$(basename \$d)".fastq"
                fi
            fi
        fi
    done
    """
}

process merge {

    publishDir "${params.out_dir}/01.basecalled/", mode: "copy"
    
    input: 
    path fq from Channel.fromPath(params.in_dir+"/basecalled").collect()

    when:
    !params.no_merge && !params.basecall

    output:
    path "merged/*.fastq" into fastq_wo_guppy

    """
    mkdir merged
    IFS=','
    read -ra BCS <<<"${params.barcodes}"
    # als maar 1 barcode meegegeven wordt met leading 0: octaal voor bash
    if [[ \${#BCS[@]} -eq 1 ]] && [[ \${#BCS[0]} -eq 1 ]]; then
        BCS[0]="0"\${BCS[0]}
    fi

    for d in ./basecalled/*
    do
        # merge without barcodes
        if [[ -f \$d ]]; then
            cat ./basecalled/*.fastq > ./merged/sample.fastq
            break
        # merge with barcodes
        elif [[ \$d =~ barcode* ]]; then
            # merge all the barcodes
            if [[ \${#BCS[@]} -eq 0 ]]; then
                cat \$d"/"* > ./merged/\$(basename \$d)".fastq"
            # merge only the barcodes from params.barcodes
            else
                bc=\${d##*barcode}
                if [[ " \${BCS[@]} " =~ " \$bc " ]]; then
                    cat \$d"/"* > ./merged/\$(basename \$d)".fastq"
                fi
            fi
        fi
    done
    """
}

process no_merge {

    publishDir "${params.out_dir}/01.basecalled/merged", mode: "copy"
    
    input: 
    path fq from Channel.fromPath(params.in_dir+"/basecalled/merged/*")

    when:
    params.no_merge && !params.basecall

    output:
    path "*.fastq" into fastq_merged

    """
    IFS=','
    read -ra BCS <<<"${params.barcodes}"
    # als maar 1 barcode meegegeven wordt met leading 0: octaal voor bash
    if [[ \${#BCS[@]} -eq 1 ]] && [[ \${#BCS[0]} -eq 1 ]]; then
        BCS[0]="0"\${BCS[0]}
    fi

    if [[ \${#BCS[@]} -ne 0 ]]; then
	    file=${fq}
	    basename=\${file%.*}
	    bc=\${basename##barcode}
	    if [[ " \${BCS[@]} " =~ " \$bc " ]]; then
            echo "\$bc"
 	        find . -type f -name "barcode\${bc}*.fastq" > barcode\${bc}*.fastq
	    fi
    fi
    """

}

fastqs = fastq_w_guppy.mix(fastq_wo_guppy, fastq_merged)
fastqs.into{fastq_nanoplot;fastq_nanocomp;fastq_assembly}

process nanoplot {

    container 'staphb/nanoplot:1.33.0'

    publishDir "${params.out_dir}/02.qc/nanoplot", mode:"copy"

    input:
    // .flatten() analyses all the barcodes separately 
    path fq from fastq_nanoplot

    output:
    path "*"

    when:
    params.nanoplot

    script:
    """
    NanoPlot -t 4 -o ./nanoplot --fastq_rich $fq
    """
}

process nanocomp {

    container 'bikc/nanocomp:v1'

    publishDir "${params.out_dir}/02.qc/nanocomp", mode:"copy"

    input:
    // .flatten() analyses all the barcodes separately
    path fq from fastq_nanocomp

    output:
    path "*"
    
    when:
    params.nanocomp 

    script:
    """
    NanoComp -t 4 -o ./nanocomp -f png --plot violin --fastq $fq
    """
}

process flye_assembly {
        
    container 'staphb/flye:2.8'
    
    errorStrategy 'ignore'
        
    // pattern: publish everything except for the fqs
    publishDir "${params.out_dir}/03.assembly/", mode: "copy"

    input: 
    path fq from fastq_assembly.flatten()

    output:
    // **: search in above directories
    path "*/*assembly*" optional true
    path "**flye.log" optional true 
    tuple file(fq), file("*/*assembly.fasta") optional true into (flye_prokka, mapping)

    when:
    params.assemble

    script:

    gsize = params.gsize ? "--genome-size $params.gsize" : ""
    plasmids = params.plasmids ? "--plasmids" : ""
    coverage = params.asm_coverage ? "--asm-coverage $params.asm_coverage" : ""
    meta = params.meta ? "--meta" : ""

    """
    flye --nano-raw $fq $gsize $coverage $plasmids $meta --out-dir ./${fq.simpleName} -t ${params.assemblyP} -i 5

    # add bc on assembly.fasta
	for file in ./${fq.simpleName}/assembly*; do
		basename=\${file##*/}
		basename_no_ext=\${basename%.*}
		extension=\${file##*.}
		mv "\$file" "./${fq.simpleName}/${fq.simpleName}_\$basename_no_ext.\$extension"
	done

    """
} 

process mapping {

    container 'staphb/minimap2:2.17'

    publishDir "${params.out_dir}/04.mapping/", mode: "copy", pattern: "*.sam"

    input:
    tuple file(fq), file(assembly) from mapping

    output:
    path "*" into sam_files
    tuple file(fq), file(assembly), file("*.sam") into racon

    when:
    params.mapping || params.polishing
    
    """
    minimap2 -a $assembly $fq > ${fq.simpleName}.sam
    """
}

process sam_to_bam {

    container 'staphb/samtools:1.11'

    publishDir "${params.out_dir}/04.mapping", mode: "copy"

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

    publishDir "${params.out_dir}/04.mapping", mode: "copy"

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

    publishDir "${params.out_dir}/04.mapping", mode: "copy"

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
    
    /* when data is limited, sometimes Flye produces an assembly empty assembly, which will cause an error.
    errorStragy 'ignore' will not terminate te script if this occurs */

    errorStrategy 'ignore'    
   
    container 'staphb/racon:1.4.20'

    publishDir "${params.out_dir}/05.polishing/racon", mode: "copy" , pattern: "*_racon.fasta"

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

    container 'bikc/medaka:v1'

    publishDir "${params.out_dir}/05.polishing/medaka", mode: "copy", pattern: "**.fasta"

    input:
    tuple file(fq), file(assembly) from medaka

    output:
    path "**.fasta" into medaka_result
    tuple file(fq), file ("**consensus.fasta") into medaka_prokka

    when:
    params.polishing

    """
    medaka_consensus -i $fq -d $assembly -o ./${fq.simpleName} -t ${params.threads_polishing} -m ${params.model} 
    mv ./${fq.simpleName}/consensus.fasta ${fq.simpleName}_consensus.fasta
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
    
    publishDir "${params.out_dir}/06.annotation/", mode: "copy"

    input:
    tuple file(fq), file(assembly) from prokka_tuple

    output:
    path "*" into prokka_result

    when:
    params.annotation

    """
    prokka --outdir ./${fq.simpleName} --prefix ${fq.simpleName} $assembly --cpus ${params.threads_annotation}
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
   	default: false
    • barcode_kits:                 | ${c_yellow}If provided, this will demultiplex the samples using the provided kit${c_reset}
    • bc_config:                    | ${c_yellow}If provided, will use a different config file then default${c_reset}
        default: dna_r9.4.1_450bps_fast.cfg
    • skip_qc:                      | ${c_yellow}If provided, this will not split basecalled reads into pass or fail${c_reset}
   	default: false 
    • num_callers (default:1):      | ${c_yellow}Number of callers to use for basecalling (1=4 threads!)${c_reset}
	default: 1
    ${c_cyan}Merging:${c_reset}  
    • no_merge:                     | ${c_yellow}If provided: will expect a /basecalled/ folder with fastq files for each barcode / a sample${c_reset}
    	default: false              | ${c_yellow}If not provided: will expect a /basecalled/merged/ folder with merged fastq files for each barcode${c_reset}
                                    | ${c_yellow}/ a sample${c_reset}
   
    ${c_cyan}QC:${c_reset}                             | ${c_blue}Quality Control steps${c_reset}
    • nanocomp:                     | ${c_yellow}If provided, will perform nanocomp analysis ${c_reset}
	default: true  
    • nanoplot:                     | ${c_yellow}If provided, will perform nanoplot analysis ${c_reset}
	default: true

    ${c_cyan}Assembly:${c_reset}                       | ${c_blue}Only when providing --assemble option${c_reset}
    • assemble:                     | ${c_yellow}If provided, this will assemble the genomes using Flye ${c_reset}
	default:true
    • gsize:                        | ${c_yellow}Expected genome size (not mandatory) ${c_reset}
    • meta:                         | ${c_yellow}Metagenome / Uneven coverage${c_reset}
    • plasmids:                     | ${c_yellow}rescue short unassembled plasmids${c_reset}
    • asm_coverage:                 | ${c_yellow}reduced coverage for initial disjointig assembl${c_reset}
    • assemblyP	 	                | ${c_yellow}Number of threads per barcode (use max: 32/nbarcodes) ${c_reset}
	default: 8   
    • barcodes:                     | ${c_yellow}Comma separated list of barcode numbers that are expected, if barcodes are provided${c_reset}
                                    | ${c_yellow}Numbers should include the leading 0s. E.g. 03,08,11${c_reset}

    ${c_cyan}mapping:${c_reset}                        | ${c_blue}If provided, will map, sort and index sequences ${c_reset}
	default: true

    ${c_cyan}polishing:${c_reset}                      | ${c_blue}If provided, will polish sequences (requires mapping) ${c_reset}
	default: true   
    • threads_polishing             | ${c_yellow}Number of threads used for polishing ${c_reset}
   	default: 8
    • model:                        | ${c_yellow}Model used for Medaka polishing: {pore}_{device}_{caller variant}_{caller version} ${c_reset}
        default: r941_min_fast_g303
	    hac: r941_min_high_g360
	    bonito: bonito_medaka_0.3.1.hdf5

    ${c_cyan}annotation:${c_reset}                     | ${c_blue}If provided, will anotate sequences (requires mapping) ${c_reset}
   	default: true
    • threads_annotation	    | ${c_yellow}Number of threads used for annotation ${c_reset} 
        default: 8  
 
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
    • meta:                         | ${c_yellow}${params.meta}${c_reset}
    • plasmids:                     | ${c_yellow}${params.plasmids}${c_reset}
    • asm_coverage:                 | ${c_yellow}${params.asm_coverage}${c_reset}
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
    • polishing:                    |"""

    if (params.polishing) {
        text = text + """ ${c_yellow}${params.polishing}${c_reset}
    • model:                        | ${c_yellow}${params.model} ${c_reset}
    • threads_polishing:            | ${c_yellow}${params.threads_polishing} ${c_reset}

    """
    }
    else {
        text = text + """ ${c_yellow}${params.polishing}${c_reset}

    """
    }

    text = text + """${c_cyan}Annotation:${c_reset}
    • annotation:                   | """

    if (params.annotation) {
	text = text + """${c_yellow}${params.annotation}${c_reset}
    • threads_annotation:           | ${c_yellow}${params.threads_annotation}${c_reset} 

    """
    }
    else {
	text = text + """${c_yellow}${params.annotation}${c_reset}
    
    """
    } 
   

    text=text+"""Output directory:  | ${c_yellow}${params.out_dir}${c_reset}
    =================================================="""


    return text;
}

