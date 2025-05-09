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

//////////////////////////////////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////
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

/* Quality control processes */
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

/* Trimming process */
process trimming{
    label "trimming"

    publishDir "${params.out_dir}/03.trimming/", mode: "copy"

    input:
    path fq

    output:
    path "${fq.simpleName}_trimmed.fastq"

    script:
    """
    pychopper -t ${params.t_trimming} ${fq} ${fq.simpleName}_trimmed.fastq
    # -k kit PCS109, PCS110, PCS111, LSK114
    """
}

/* Mapping to human reference genome */
process human_mapping {

    label "mapping"

    publishDir "${params.out_dir}/04.mapped_reads/human", mode: "copy"

    input:
    path fq             // Input FASTQ files either from trimming or preprocessing

    output:
    path "*_human.sam", emit: sam_output
    path "*_unmapped_human.sam", emit: unmapped_output

    script:
    """
    minimap2 -ax splice --secondary=no ${params.human_ref} $fq -t ${params.t_mapping} > ${fq.simpleName}_human.sam
    """
}

/* Extract unmapped/mapped reads and sort/index the mapped reads */
process sort_index_human {

    label "samtools"

    publishDir "${params.out_dir}/04.mapped_reads/sorted_indexed", mode: "copy"

    input:
    path mapped_human       // Mapped sequences from human genome

    output:
    path "*.sorted.bam"
    path "*.sorted.bam.bai"

    script:
    """
    # Extract unmapped reads
    samtools view -h -f 4 ${fq.simpleName}_human.sam > ${fq.simpleName}_unmapped_human.sam
    
    # Extract mapped reads
    samtools view -h -F 4 ${fq.simpleName}_human.sam > ${fq.simpleName}_mapped_human.sam

    # Sort and index BAM files
    samtools sort -o ${mapped_human.simpleName}_mapped_sorted.bam ${mapped_human}
    samtools index ${mapped_human.simpleName}_mapped_sorted.bam
    """
}

/* Mapping of unmapped human reads to candida reference genome */
process candida_mapping {

    label "mapping"

    publishDir "${params.out_dir}/04.mapped_reads/candida", mode: "copy"

    input:
    path unmapped_human     // Input FASTQ file from human_mapping process

    output:
    path "*_candida.sam", emit: sam_output
    path "*_unmapped_candida.sam", emit: unmapped_output

    script:
    """
    minimap2 -ax splice --secondary=no ${params.candida_ref} $unmapped_human -t ${params.t_mapping} > ${unmapped_human.simpleName}_candida.sam
    """
}

/* Extract unmapped/mapped reads and sort/index the mapped reads */
process sort_index_candida {

    label "samtools"

    publishDir "${params.out_dir}/04.mapped_reads/sorted_indexed", mode: "copy"

    input:
    path mapped_candida     // Mapped sequences from candida genome

    output:
    path "*.sorted.bam"
    path "*.sorted.bam.bai"

    script:
    """
    # Extract unmapped reads
    samtools view -h -f 4 ${fq.simpleName}_candida.sam > ${fq.simpleName}_unmapped_candida.sam
    
    # Extract mapped reads
    samtools view -h -F 4 ${fq.simpleName}_candida.sam > ${fq.simpleName}_mapped_candida.sam

    # Sort and index BAM files
    samtools sort -o ${mapped_candida.simpleName}_mapped_sorted.bam ${mapped_candida}
    samtools index ${mapped_candida.simpleName}_mapped_sorted.bam
    """
}

/* Mapping of unmapped candida reads to bacterial reference genome */
process bacterial_mapping {

    label "mapping"

    publishDir "${params.out_dir}/04.mapped_reads/bacteria", mode: "copy"

    input:
    path unmapped_candida       // unmapped candida reads

    output:
    path "*_bacterial.sam", emit: sam_output
    path "*_unmapped_bacteria.sam", emit: unmapped_output

    script:
    """
    minimap2 -ax map-ont --secondary=no ${params.bacteria_index} $unmapped_candida -t ${params.t_mapping} > ${unmapped_candida.simpleName}_bacteria.sam
    """
}

/* Sort and index the mapped reads */
process sort_index_bacteria {

    label "samtools"

    publishDir "${params.out_dir}/04.mapped_reads/sorted_indexed", mode: "copy"

    input:
    path mapped_bacteria    // Mapped sequences from bacteria index

    output:
    path "*.sorted.bam"
    path "*.sorted.bam.bai"

    script:
    """
    # Extracting unmapped reads & transforming it to fastq for Kraken
    samtools view -b -f 4 ${mapped_bacteria.simpleName}_bacteria.sam | samtools fastq - > ${mapped_bacteria.simpleName}_unmapped_bacteria.fastq

    # Extract mapped reads
    samtools view -h -F 4 ${fq.simpleName}_bacterial.sam > ${fq.simpleName}_mapped_bacteria.sam

    # Sort and index BAM files
    samtools sort -o ${mapped_bacteria.simpleName}_mapped_sorted.bam ${mapped_bacteria}
    samtools index ${mapped_bacteria.simpleName}_mapped_sorted.bam
    """
}

/* Kraken implementation for left over unmapped reads*/
process kraken_classification {

    label "kraken"

    publishDir "${params.out_dir}/05.kraken_results/", mode: "copy"

    input:
    path unmapped_bacterial_fastq,  // Input: Unmapped bacterial reads
    path kraken_db                  // Kraken database

    output:
    path "*.kraken", emit: kraken_output

    script:
    """
    kraken2 --db ${kraken_db} --report ${unmapped_bacterial_fastq.simpleName}_kraken_report.txt --output ${unmapped_bacterial_fastq.simpleName}_kraken_output.txt ${unmapped_bacterial_fastq}
    """
}

/* Count genes with NanoCount */
process quantification {
    
    label "counting"

    publishDir "${params.out_dir}/06.quantification/", mode: "copy"

    input:
    path mapped_human
    path mapped_candida
    path mapped_bacteria

    output:
    path "*.counts"

    script:
    """
    # Counting human genes
    featureCounts -T ${params.t_counting} -a ${params.human_annotation} -o human_gene_counts.txt --primary ${mapped_human}

    # Counting candida genes
    featureCounts -T ${params.t_counting} -a ${params.candida_annotation} -o candida_gene_counts.txt --primary ${mapped_candida}

    # Counting bacterial genes
    featureCounts -T ${params.t_counting} -a ${params.bacteria_annotation} -o bacterial_gene_counts.txt --primary ${mapped_bacteria}

    # flags
    # -T: number of threads
    # -a: annotation file
    # -o: output file
    # -t: feature type in GTF file (default: exon, gene, CDS, UTR or Transcripts)
    # -g: attribute type (default: gene_id, transcript_id, exon_id, gene_name, biotype)
    # --primary: only primary alignments
    """
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Workflow //
workflow {
    /* Create channel to start from */
    fastq_files = channel
        .fromPath("${params.in_dir}/fastq/")
        .collect()

    preprocessed_fastq = preprocess(fastq_files)

    /// QUALITY CONTROL ///
    if (params.qc) {
        nanoplot_QC(preprocessed_fastq.flatten())
        nanocomp_QC(preprocessed_fastq.collect())
    }

    // Trimming
    if (params.trimming) {
        trimmed_fastq = trimming(preprocessed_fastq.flatten())
        mapping_input = trimmed_fastq
    } else {
        mapping_input = preprocessed_fastq.flatten()
    }

    // Mapping to human reference genome
    human_sam_files = human_mapping(mapping_input)

    // Mapping to Candida reference genome
    candida_sam_files = candida_mapping(human_sam_files.unmapped_output)

    // Mapping to bacterial reference genome
    bacterial_sam_files = bacterial_mapping(candida_sam_files.unmapped_output)

    // Kraken classification (using unmapped bacterial reads in FASTQ format)
    kraken_classification(
        sort_index_bacteria.unmapped_bacteria_fastq,
        params.kraken_db
    )

    // Quantification
    quantification(
        mapped_human.sam_output,
        mapped_candida.sam_output,
        mapped_bacteria.sam_output,
        params.human_annotation,
        params.candida_annotation,
        params.bacteria_annotation
    )
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
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
     ====================================================================================================
     Provided parameters:
 
     ${c_cyan}General:${c_reset} 
     • in_dir:                                  | ${c_yellow}${params.in_dir}${c_reset}
     • out_dir:                                 | ${c_yellow}${params.out_dir}${c_reset}
     • barcodes:                                | ${c_yellow}${params.barcodes} ${c_reset}
     • ignore_errors                            | ${c_yellow}${params.ignore_errors}${c_reset}
     • help:                                    | ${c_yellow}${params.help}${c_reset}
     • monochrome_logs:                         | ${c_yellow}${params.monochrome_logs}${c_reset}

     ${c_cyan}QC:${c_reset}
     • qc:                                      | ${c_yellow}${params.qc}${c_reset}
     • t_qc:                                    | ${c_yellow}${params.t_qc} ${c_reset}

     ${c_cyan}trimming:${c_reset}
     • trimming                                 | ${c_yellow}${params.trimming}${c_reset}
     • t_trimming                               | ${c_yellow}${params.t_trimming}${c_reset}                

     ${c_cyan}mapping:${c_reset}                | ${c_blue}${params.mapping}${c_reset}
     • t_mapping                                | ${c_yellow}${params.t_mapping}${c_reset}

     ${c_cyan}reference genomes:${c_reset}      | ${c_blue}${params.mapping}${c_reset}
     • human_ref                                | ${c_yellow}${params.human_ref}${c_reset}
     • candida_ref                              | ${c_yellow}${params.candida_ref}${c_reset}
     • bacteria_index                           | ${c_yellow}${params.bacteria_index}${c_reset}

    ${c_cyan}Kraken classification:${c_reset}   | ${c_blue}${params.mapping}${c_reset}

     ${c_cyan}Quantification:${c_reset}         | ${c_blue}${params.mapping}${c_reset}
     ====================================================================================================
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
    ====================================================================================================
    Parameters:

    ${c_cyan}General:${c_reset}                         | ${c_red}Must supply parameters${c_reset}
    • in_dir:                                           | ${c_yellow}Input directory (directory with fastq folder)${c_reset}
    • out_dir:                                          | ${c_yellow}Output directory for results${c_reset}
                                                        | ${c_red}Optional parameters${c_reset}
    • barcodes:                                         | ${c_yellow}Comma separated list of barcode numbers that the user wants to analyse if barcodes are present${c_reset}
                                                        | ${c_yellow}All barcodes are automatically analysed if barcodes are present and the --barcodes parameter is not provided${c_reset}
                                                        | ${c_yellow}Numbers should include the leading 0s. E.g. 03,08,11${c_reset}
    • ignore_errors                                     | ${c_yellow}Option to continue running the pipeline even if errors occur ${c_reset} 
    default: false
    • help                                              | ${c_yellow}Show this${c_reset}
    default: false   
    • monochrome_logs                                   | ${c_yellow}If set to true, the logs will be printed in monochrome${c_reset}
    default: false

    ${c_cyan}QC:${c_reset}                              | ${c_blue}Quality Control steps${c_reset}
    • qc:                                               | ${c_yellow}If provided, will perform QC analysis (NanComp and NanoPlot) ${c_reset}
	default: true  
    • t_qc                                              | ${c_yellow}Number of threads used for QC${c_reset}
	default: 4

    ${c_cyan}trimming:${c_reset}                        | ${c_blue}Trimming steps${c_reset}
    • trimming:                                         | ${c_yellow}If provided, will perform trimming (Porechop) ${c_reset}
	default: false
    • t_trimming                                        | ${c_yellow}Number of threads used for trimming${c_reset}
	default: 4

    ${c_cyan}mapping:${c_reset}                         | ${c_blue}If provided, will map, sort and index sequences ${c_reset}
	default: true
    • t_mapping                                         | ${c_yellow}Number of threads used for mapping ${c_reset}
    default: 4

    ${c_cyan}reference genomes:${c_reset}               | ${c_blue}Reference genomes${c_reset}
    • human_ref                                         | ${c_yellow}Path to the human reference genome${c_reset}
    • candida_ref                                       | ${c_yellow}Path to the candida reference genome${c_reset}
    • bacteria_index                                    | ${c_yellow}Path to the bacteria index${c_reset}

    ${c_cyan}Kraken classification:${c_reset}           | ${c_blue}If provided, will perform quantification ${c_reset}
    default: true 

    ${c_cyan}Quantification:${c_reset}                  | ${c_blue}If provided, will perform quantification ${c_reset}
    default: true   

    ====================================================================================================                    
    Run using: nextflow run main.nf --in_dir {Input directory} --out_dir {Output directory}
    ====================================================================================================
    """
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
