#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Help message //
if (params.help) {
    print helpMessage()
    exit 0
}

else {
    print parameterShow()
}

// Parameters validation //
/* Check if there is an input directory specified */
error = false
if (!params.in_dir){
    println "${c_red}\nParameter --in_dir is missing${c_reset}"
    error = true
} 

/* Check if input directory has files */
else if (params.in_dir instanceof Boolean) {
    println "${c_red}\nParameter --in_dir is empty, no input directory specified${c_reset}"
    error = true
} 

/* Check if input directory exists */
else if (!file(params.in_dir).exists()){
    println "${c_red}\nParameter --in_dir is not an existing directory${c_reset}"
    error = true
}

/* Check if there is an output directory specified */
if (!params.out_dir){
    println "${c_red}\nParameter --out_dir is missing${c_reset}"
    error = true
} 

/* Check if output directory has files */
else if (params.out_dir instanceof Boolean) {
    println "${c_red}\nParameter --out_dir is empty, no output directory specified${c_reset}"
    error = true
} 

/* Check if output directory exists */
else if (!file(params.out_dir).isDirectory()){
    println "${c_red}\nParameter --out_dir is not an existing directory${c_reset}"
    error = true
}

/* Check if there is a fastq directory inside the input directory that contains the fastq files */
fastqDir = file("${params.in_dir}/fastq/")
if (!fastqDir.isDirectory()) {
    println "${c_red}\nParameter --in_dir needs to contain a fastq directory${c_reset}"
    println "${c_red}Make a fastq directory inside the in_dir${c_reset}"
    error=true
} 

/* Check if there only subdirectories in the fastq directory, no files */
else {
    containsFiles = false
    fastqDir.eachFile { file ->
        if (file.isFile()) {
            containsFiles = true
        }
    }

    if (containsFiles) {
        println "${c_red}\nThe fastq directory should only contain subdirectories (with barcodes or 1 sample), no files are allowed.${c_reset}"
        error = true
    }
}

/* Check if the fastq input directory contains fastq files or fastq.gz files*/
if (file("${params.in_dir}/fastq/").isDirectory()) {
    count = 0
    file("${params.in_dir}/fastq").eachFileRecurse {
        if (it.name.endsWith(".fastq") || it.name.endsWith(".fastq.gz")) {
            count += 1
        }
    }
    if (count == 0) {
        println "${c_red}\nNo fastq or fastq.gz files found in the fastq directory${c_reset}"
        error=true
    }
}

/* Check if the barcodes that the user wants to analyze, are found in the subdirectories of the fastq directory */
if (params.barcodes) {
    missingBarcodes = [] 
    params.barcodes.each { barcode ->
        String formattedBarcode = String.format("%02d", barcode.toInteger())
        barcodeDir = file("${params.in_dir}/fastq/barcode${formattedBarcode}")
        if (!barcodeDir.isDirectory()) {
            missingBarcodes << formattedBarcode 
        }
    }

    if (missingBarcodes) {
        println "${c_red}\nFollowing barcodes are not found in the 'in_dir/fastq/' directory: ${missingBarcodes.join(', ')}${c_reset}"
        error = true
    }
}

if (error) {
    exit 0
}


// processes //
/* preprocess data involves merging the fastq files per barcode and unzipping the FASTQ files if necessary */
process preprocess {

    publishDir "${params.out_dir}/01.preprocessed_data/", mode: "copy"

    input:
    path fastq_files
    
    output:
    path '*.fastq' 

    script:
    """
    # Read barcodes from parameters
    IFS=',' read -ra barcodes <<<"${params.barcodes}"

    # Handle barcodes to add leading zeros if necessary (f. ex --barcodes 4)
    for i in "\${barcodes[@]}"; do
        if [[ \${#i} -eq 1 ]]; then
            i="0\$i"  # Add leading zero
        fi
        leading_barcodes+=("\$i")  # Add to new array
    done

    # Check for barcodes in directories
    for dir in ./fastq/*; do
        # Check if the subdirectory is a directory (f. ex barcode13)
        if [[ -d "\$dir" ]]; then
            dir_name=\$(basename "\$dir")
            # Check if the run was barcoded or not
            if [[ "\$dir" == *"barcode"* ]]; then
                 # Check if the user want to analyze specific barcodes
                 # If the user did not specify any barcodes, all the barcode directories are analyzed
                if [[ -z "${params.barcodes}" ]]; then
                    # Check if there are only fastq.gz files in the directory
                    if ls "\$dir"/*.fastq.gz 1> /dev/null 2>&1; then
                        cat \$dir"/"* | gunzip -c > ./\$dir_name".fastq"
                    # Check if there are only fastq files in the directory
                    elif ls "\$dir"/*.fastq 1> /dev/null 2>&1; then
                        cat \$dir"/"*  > ./\$(dir_name)".fastq"
        
                    fi
                else
                    # Find the barcodes that the user wants to analyze
                    bc=\${dir_name##*barcode}
                    if [[ "\${leading_barcodes[@]}" =~ "\$bc" ]]; then
                        # Check if there are fastq.gz files in the direcotry
                        if ls "\$dir"/*.fastq.gz 1> /dev/null 2>&1; then
                            cat \$dir"/"* | gunzip -c > ./\$dir_name".fastq"
                        # Check if there are fastq files in the directory
                        elif ls "\$dir"/*.fastq.gz 1> /dev/null 2>&1; then
                            cat \$dir"/"*  > ./\$dir_name".fastq"
                        fi
                    fi
                fi
            else
                # If there are no barcodes, then the directory basename is used as name
                # Check if there are fastq.gz files in the directory
                if ls "\$dir"/*.fastq.gz 1> /dev/null 2>&1; then
                    cat \$dir"/"* | gunzip -c > ./\$dir_name".fastq"
                # Check if there are fastq files in the directory
                elif ls \$dir"/"* 1> /dev/null 2>&1; then
                    cat \$dir"/"* > ./\$dir_name".fastq"
                fi
            fi
        fi
    done
    """
}

process nanoplot_QC {

    label "qc"

    publishDir "${params.out_dir}/02.qc/nanoplot/${fq.simpleName}", mode:"copy", pattern: "*.{png,txt,html}"

    input:
    path fq

    output:
    path "*" 

    script:
    """
    NanoPlot -t ${params.t_qc} -o . --N50 --fastq_rich ${fq}
    """
}

process nanocomp_QC {

    label "qc"

    publishDir "${params.out_dir}/02.qc/nanocomp/", mode:"copy"

    input:
    path fq

    output:
    path "*" 

    script:
    """
    NanoComp -t ${params.t_qc} -o . -f png --plot violin --fastq ${fq}
    """
}

//process flye_assembly {
        
//    label "flye"
        
//    publishDir "${params.out_dir}/03.flye_assembly/", mode: "copy", pattern: "**.{fasta,gfa,gv,txt,log}"

//    input: 
//    path fq

//    output:
//    path "**.log" 
//    path "*/*assembly*" 
//    tuple file(fq), file("*/*assembly.fasta"), emit: fasta_flye
//    path("*/*assembly.fasta"), emit: fasta_flye_annot

//    script:

//    gsize = params.gsize ? "--genome-size $params.gsize" : ""
//    coverage = params.asm_coverage ? "--asm-coverage $params.asm_coverage" : ""
//    meta = params.meta ? "--meta" : ""
//    read_quality = params.nano_hq ? "--nano-hq" : "--nano-raw"

//    """
//    flye $read_quality $fq $gsize $coverage $meta --out-dir ./${fq.simpleName} -t ${params.t_assembly} -i 5

//    # add bc on assembly.fasta
//	for file in ./${fq.simpleName}/assembly*; do
//		basename=\${file##*/}
//		basename_no_ext=\${basename%.*}
//		extension=\${file##*.}
//		mv "\$file" "./${fq.simpleName}/${fq.simpleName}_\$basename_no_ext.\$extension"
//	done
//    """
//}

process initial_alignment {

    label "mapping"

    publishDir "${params.out_dir}/04.mapped_reads/inital", mode: "copy"

    input:
    tuple file(fq), file(assembly)

    output:
    path "*.sam", emit: sam_output

    """
    minimap2 -ax map-ont $assembly $fq -t ${params.t_mapping} > ${fq.simpleName}.sam
    """
}

process initial_bam_processing {

    label "samtools"

    publishDir "${params.out_dir}/04.mapped_reads/inital", mode: "copy", pattern: "*.{bam,bai}"

    input:
    path sam_file 

    output:
    path "*"

    """
    # Convert SAM to BAM, sort and index
    samtools view -S -b $sam_file | samtools sort -@ ${params.t_mapping} -o ${sam_file.simpleName}_sorted.bam
    samtools index ${sam_file.simpleName}_sorted.bam
    """
}


//process medaka_polishing {

//    label "medaka"

//    publishDir "${params.out_dir}/05.polishing/medaka", mode: "copy", pattern: "**.fasta"

//    input:
//    tuple file(fq), file(assembly) 

//    output:
//    tuple file(fq), file ("**consensus.fasta"), emit: fasta_medaka
//    path ("**consensus.fasta"), emit: fasta_medaka_annot

//    script:

//    model = params.model ? "-m $params.model" : ""

//    """
//    medaka_consensus -i $fq -d $assembly -o ./${fq.simpleName} -t ${params.t_polishing} $model 
//    mv ./${fq.simpleName}/consensus.fasta ${fq.simpleName}_consensus.fasta
//    """

//}

process remapped_alignment {

    label "mapping"

    publishDir "${params.out_dir}/04.mapped_reads/remapped", mode: "copy"

    input:
    tuple file(fq), file(assembly)

    output:
    path "*.sam", emit: sam_output

    """
    minimap2 -ax map-ont $assembly $fq -t ${params.t_mapping} > ${fq.simpleName}.sam
    """
}

process remapped_bam_processing {

    label "samtools"

    publishDir "${params.out_dir}/04.mapped_reads/remapped", mode: "copy", pattern: "*.{fasta,bam,bai}"

    input:
    path sam_file 

    output:
    path "*"

    """
    # Convert SAM to BAM, sort and index
    samtools view -S -b $sam_file | samtools sort -@ ${params.t_mapping} -o ${sam_file.simpleName}_sorted.bam
    samtools index ${sam_file.simpleName}_sorted.bam
    """
}


//process prokka_annotation {
  
//    label "prokka"
    
//    publishDir "${params.out_dir}/06.annotation/prokka", mode: "copy"

//   input:
//    path(assembly) 

//   output:
//    path "*"

//    """
//    prokka --outdir ./${assembly.simpleName} --prefix ${assembly.simpleName} $assembly --cpus ${params.t_fast_annotation}
//    """
//}

//process pgap_annotation {

//    label "pgap"

//    publishDir "${params.out_dir}/06.annotation/pgap", mode: "copy"

//    input:
//    path(assembly)

//    output:
//    path "*"

//    script:
//    errors = params.ignore_errors ? "--ignore-all-errors" : ""

//    """
//    cp ${assembly} ${params.organism}.fasta

//    create_yaml.py ${params.organism}

//    python ${params.pgap_path}/pgap.py -n $errors --no-self-update -o ./PGAP --use-version ${params.pgap_version} ./${params.organism}.yaml --cpus ${params.t_accurate_annotation}

//    """

//}


//process busco_qc{

//    label "busco"
//    containerOptions '-u $(id -u):$(id -g)'
        
//    publishDir "${params.out_dir}/07.genome_qc/", mode: "copy"

//    input:
//    tuple file(fq), file(assembly), path(data) 

//    output:
//    path "*"

//    """
//    busco -m genome -i $assembly -o ./${fq.simpleName} -l ${params.lineage} --offline --download_path ${data} --cpu ${params.t_assembly_qc}
//    """

//}




// Workflow //
workflow {
    /* Create channel to start from */
    fastq_files = channel
        .fromPath("${params.in_dir}/fastq/")
        .collect()

    preprocessed_fastq = preprocess(fastq_files)

    if (params.qc) {
        nanoplot_QC(preprocessed_fastq.flatten())
        nanocomp_QC(preprocessed_fastq.collect())
    }

    /*draft_assembly=flye_assembly(preprocessed_fastq.flatten())

    if (params.polishing) {
        polished_assembly = medaka_polishing(draft_assembly.fasta_flye)
        if (params.mapping) {
            initial_sam_files = initial_alignment(draft_assembly.fasta_flye)
            initial_bam_files = initial_bam_processing(initial_sam_files.sam_output)
            remapped_sam_files = remapped_alignment(polished_assembly.fasta_medaka)
            remapped_bam_files = remapped_bam_processing(remapped_sam_files.sam_output)
        }
        if (params.fast_annotation) {
            prokka = prokka_annotation(polished_assembly.fasta_medaka_annot)
        }

        if (params.accurate_annotation) {
            pgap = pgap_annotation(polished_assembly.fasta_medaka_annot)
        }

        if (params.assembly_qc) {
            busco_tuple = polished_assembly.fasta_medaka.combine(Channel.fromPath(params.busco_path))
            busco = busco_qc(busco_tuple)
        }
    }
    if (!params.polishing) {
        if (params.mapping) {
            initial_sam_files = initial_alignment(draft_assembly.fasta_flye)
            initial_bam_files = initial_bam_processing(initial_sam_files.sam_output)
        }
        if (params.fast_annotation) {
            prokka = prokka_annotation(draft_assembly.fasta_flye_annot)
        }

        if (params.accurate_annotation) {
            pgap = pgap_annotation(draft_assembly.fasta_flye_annot)
        }

        if (params.assembly_qc) {
            busco_tuple = draft_assembly.fasta_flye.combine(Channel.fromPath(params.busco_path))
            busco = busco_qc(busco_tuple)
        } 
    } */
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
    c_gray = params.monochrome_logs ? '' : "\033[0;241m"


    return """\
     ==================================================
     Provided parameters:
 
     ${c_cyan}General:${c_reset} 
     • in_dir:                       | ${c_yellow}${params.in_dir}${c_reset}
     • out_dir:                      | ${c_yellow}${params.out_dir}${c_reset}
     • barcodes:                     | ${c_yellow}${params.barcodes} ${c_reset}


     ${c_cyan}QC:${c_reset}
     • qc:                           | ${c_yellow}${params.qc}${c_reset}
     • t_qc:                         | ${c_yellow}${params.t_qc} ${c_reset}


    ${c_cyan}assembly:${c_reset}                      
    • nano_hq:                      | ${c_yellow}${params.nano_hq}${c_reset}
    • gsize:                        | ${c_yellow}${params.gsize}${c_reset}
    • meta:                         | ${c_yellow}${params.meta}${c_reset}
    • asm_coverage:                 | ${c_yellow}${params.asm_coverage}${c_reset}
    • t_assembly                    | ${c_yellow}${params.t_assembly}${c_reset}

    ${c_cyan}mapping:${c_reset}                         | ${c_blue}${params.mapping}${c_reset}
    • t_mapping                     | ${c_yellow}${params.t_mapping}${c_reset}


    ${c_cyan}polishing:${c_reset}                       | ${c_blue}${params.polishing}${c_reset}
    • t_polishing                   | ${c_yellow}${params.t_polishing}${c_reset}
    • model:                        | ${c_yellow}${params.t_polishing}${c_reset}

    ${c_cyan}annotation:${c_reset}    
    • fast_annotation               | ${c_yellow}${params.fast_annotation}${c_reset}
    • t_fast_annotation             | ${c_yellow}${params.t_fast_annotation}${c_reset} 
    • accurate_annotation           | ${c_yellow}${params.accurate_annotation}${c_reset}
    • t_accurate_annotation         | ${c_yellow}${params.t_accurate_annotation}${c_reset} 
    • pgap_path                     | ${c_yellow}${params.pgap_path}${c_reset} 
    • pgap_version                  | ${c_yellow}${params.pgap_version}${c_reset} 
    • organism                      | ${c_yellow}${params.organism}${c_reset} 
    • ignore_errors                 | ${c_yellow}${params.ignore_errors}${c_reset} 

    ${c_cyan}assembly_qc:${c_reset}                     | ${c_blue}${params.assembly_qc}${c_reset}
    • t_assembly_qc                 | ${c_yellow}${params.t_assembly_qc}${c_reset} 
    • lineage                       | ${c_yellow}${params.lineage}${c_reset} 
    • busco_path                    | ${c_yellow}${params.busco_path}${c_reset} 

    ==================================================
     """
}
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

    ${c_cyan}General:${c_reset}                        | ${c_red}Must supply parameters${c_reset}
    • in_dir:                       | ${c_yellow}Input directory (directory with fastq folder)${c_reset}
    • out_dir:                      | ${c_yellow}Output directory for results${c_reset}
                                    | ${c_red}Optional parameters${c_reset}
    • barcodes:                     | ${c_yellow}Comma separated list of barcode numbers that the user wants to analyse if barcodes are present${c_reset}
                                    | ${c_yellow}All barcodes are automatically analysed if barcodes are present and the --barcodes parameter is not provided${c_reset}
                                    | ${c_yellow}Numbers should include the leading 0s. E.g. 03,08,11${c_reset}

   
    ${c_cyan}QC:${c_reset}                             | ${c_blue}Quality Control steps${c_reset}
    • qc:                           | ${c_yellow}If provided, will perform QC analysis (NanComp and NanoPlot) ${c_reset}
	default: true  
    • t_qc                          | ${c_yellow}Number of threads used for QC${c_reset}
	default: 4

    ${c_cyan}assembly:${c_reset}                       | ${c_blue}Asssembly option${c_reset}
    • nano_hq:                      | ${c_yellow}Mode for ONT Guppy5+ (SUP mode) and Q20 reads (3-5% error rate)${c_reset}
    default: true  
    • gsize:                        | ${c_yellow}Expected genome size (not mandatory) ${c_reset}
    • meta:                         | ${c_yellow}Metagenome / Uneven coverage ${c_reset}
    • asm_coverage:                 | ${c_yellow}Reduced coverage for initial disjointig assembly${c_reset}
    • t_assembly                    | ${c_yellow}Number of threads per barcode (use max: 32/nbarcodes) ${c_reset}
	default: 4   

    ${c_cyan}mapping:${c_reset}                        | ${c_blue}If provided, will map, sort and index sequences ${c_reset}
	default: true
    • t_mapping                     | ${c_yellow}Number of threads used for mapping ${c_reset}
    default: 4   

    ${c_cyan}polishing:${c_reset}                      | ${c_blue}If provided, will polish sequences (requires mapping) ${c_reset}
	default: true   
    • t_polishing                   | ${c_yellow}Number of threads used for polishing ${c_reset}
   	default: 4
    • model:                        | ${c_yellow}Model used for Medaka polishing: {pore}_{device}_{caller variant}_{caller version} f.ex. r104_e81_sup_g5015 ${c_reset}
                                    | ${c_yellow}Normally automatically selected by Medaka ${c_reset}
	    

    ${c_cyan}annotation:${c_reset}                     | ${c_blue}Annotation options${c_reset}
    • fast_annotation               | ${c_yellow}If true, Prokka will be executed ${c_reset}
    default: true
    • t_fast_annotation             | ${c_yellow}Number of threads used for Prokka annotation ${c_reset} 
    default: 4 
    • accurate_annotation           | ${c_yellow}If true, PGAP will be executed ${c_reset}
    default: false
    • t_accurate_annotation         | ${c_yellow}Number of threads used for PGAP annotation ${c_reset} 
    default: 4 
    • pgap_path                     | ${c_yellow}Path to PGAP installation files ${c_reset} 
    default: /mnt/drive3/tools/pgap/2024-07-18
    • pgap_version                  | ${c_yellow}Version of PGAP to use ${c_reset} 
    default: 2024-07-18.build7555
    • organism                      | ${c_yellow}Scientific name of the organism (necesarry for PGAP annotation) ${c_reset} 
    f.ex: organism = "escherichia_coli"
    • ignore_errors                 | ${c_yellow}Option to continue running the pipeline even if errors occur ${c_reset} 
    default: false

    ${c_cyan}assembly_qc:${c_reset}                    | ${c_blue}If provided, will calculate BUSCO scores for the assembly ${c_reset}
    default: true
    • t_assembly_qc                 | ${c_yellow}Number of threads used for BUSCO ${c_reset} 
    default: 4
    • lineage                       | ${c_yellow}Lineage to use to run BUSCO ${c_reset} 
    default: bacteria_odb10"
    • busco_path                    | ${c_yellow}Path to BUSCO installation files ${c_reset} 
    default: /data/databases/busco
 
    • help                          | ${c_yellow}Show this${c_reset}
                        
    
    Run using: nextflow run assembly.nf --in_dir {Input directory} --out_dir {Output directory}
    ==================================================
    """

}

