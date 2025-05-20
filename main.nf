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
    """
}

/* Mapping */
process mapping {
    label "mapping"

    publishDir "${params.out_dir}/04.mapping/", mode: "copy"

    input:
    path fq
    val species_list from params.MAPPING_SPECIES_LIST
    val ref_indexes from params.MAPPING_REF_INDEXES
    val mapping_params from params.MAPPING_PARAMS

    output:
    path "${fq.simpleName}_mapped_*.bam", emit: mapped_bam
    path "${fq.simpleName}_unmapped.fastq", emit: unmapped_fastq

    script:
    """
    # Start with the original FASTQ file as the input for the first mapping
    previous_unmapped="${fq}"

    # Loop through each species
    for i in "\${!species_list[@]}"; do
        species="\${species_list[\$i]}"
        ref_index="\${ref_indexes[\$i]}"
        params="\${mapping_params[\$i]}"

        echo "=== Mapping to \${species^} reference genome ==="

        # Check if the reference genome index exists
        if [ ! -f "\$ref_index" ]; then
            echo "Error: \${species^} reference genome file not found: \$ref_index"
            exit 1
        else
            echo "\${species^} reference genome file found: \$ref_index"
        fi

        # Map reads with minimap2
        sam_output="${fq.simpleName}_\${species}.sam"
        cmd="minimap2 \$params \$ref_index \$previous_unmapped > \$sam_output"
        echo "Mapping command: \$cmd"
        eval "\$cmd"
        echo "Mapped \$previous_unmapped to \${species^} genome: \$sam_output"

        # Extract unmapped reads using samtools
        unmapped_output="${fq.simpleName}_\${species}_unmapped.fastq"
        cmd_unmapped="samtools fastq -f 4 \$sam_output > \$unmapped_output"
        echo "Extraction command: \$cmd_unmapped"
        eval "\$cmd_unmapped"
        echo "Unmapped reads saved as: \$unmapped_output"

        # Sort and index the mapped reads using samtools
        sorted_bam="${fq.simpleName}_mapped_\${species}_sorted.bam"
        cmd_sort="samtools sort \$sam_output -o \$sorted_bam"
        echo "Sort command: \$cmd_sort"
        eval "\$cmd_sort"

        cmd_index="samtools index \$sorted_bam"
        echo "Index command: \$cmd_index"
        eval "\$cmd_index"
        echo "Sorted and indexed BAM: \$sorted_bam"

        # Update the input for the next iteration to be the unmapped reads from current mapping
        previous_unmapped="\$unmapped_output"
    done
    """
}

/* Kraken classification process */
process kraken_classification {
    label "kraken"

    publishDir "${params.out_dir}/05.kraken_results/", mode: "copy"

    input:
    path unmapped_fastq
    val kraken_db from params.kraken_db
    val kraken_threads from params.t_kraken

    output:
    path "*.txt", emit: kraken_reports
    path "*.fastq", emit: kraken_classified_fastq

    script:
    """
    SAMPLE_NAME=\$(basename ${unmapped_fastq} .fastq)
    SAMPLE_OUTPUT_DIR="${params.out_dir}/05.kraken_results/\${SAMPLE_NAME}"

    mkdir -p "\$SAMPLE_OUTPUT_DIR"

    kraken2 --db ${kraken_db} \
        --threads ${kraken_threads} \
        --report "\$SAMPLE_OUTPUT_DIR/\${SAMPLE_NAME}_report.txt" \
        --output "\$SAMPLE_OUTPUT_DIR/\${SAMPLE_NAME}_output.txt" \
        --use-names \
        --classified-out "\$SAMPLE_OUTPUT_DIR/\${SAMPLE_NAME}_classified.fastq" \
        --unclassified-out "\$SAMPLE_OUTPUT_DIR/\${SAMPLE_NAME}_unclassified.fastq" \
        ${unmapped_fastq}
    """
}

/* Counting process */
process counting {
    label "counting"

    publishDir "${params.out_dir}/06.quantification/", mode: "copy"

    input:
    path sample_dir
    val sample_name
    val species_list from params.COUNTING_SPECIES_LIST
    val annotations from params.COUNTING_ANNOTATIONS

    output:
    path "${sample_name}_*_bambu_results"

    script:
    """
    echo "Processing sample: ${sample_name}"

    # Loop through each species by index
    for i in "\${!species_list[@]}"; do
        SPECIES="\${species_list[\$i]}"
        ANNOTATION_FILE="\${annotations[\$i]}"

        # Find the BAM file for the species in the sample directory
        BAM_FILE="${sample_dir}/${sample_name}_mapped_\${SPECIES}_sorted.bam"
        if [ ! -f "\$BAM_FILE" ]; then
            echo "Warning: BAM file not found for species \$SPECIES in sample \$sample_name"
            exit
        fi

        OUT_DIR="${sample_name}_\${SPECIES}_bambu_results"
        mkdir -p "\$OUT_DIR"

        # Run the Bambu Runner Docker container
        docker run --rm \
            -v "${sample_dir}:/data" \
            -v "${params.out_dir}/06.quantification:/output" \
            -v "\$(dirname \$ANNOTATION_FILE):/annotations" \
            mathiasverbeke/bambu_runner:latest \
            run_bambu.R \
            --reads "/data/\$(basename \$BAM_FILE)" \
            --annotations "/annotations/\$(basename \$ANNOTATION_FILE)" \
            --genome "/annotations/genome.fa" \
            --output-dir "/output/\${sample_name}/\${SPECIES}_bambu_results" \
            --ncore ${params.t_counting} \
            --stranded no \
            --quant yes \
            --discovery no \
            --verbose yes \
            --rc-out-dir "/output/\${sample_name}/\${SPECIES}_bambu_results/rc_cache" \
            --low-memory no
    done
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

    // Mapping
    mapped_results = mapping(
        mapping_input,
        params.MAPPING_SPECIES_LIST,
        params.MAPPING_REF_INDEXES,
        params.MAPPING_PARAMS
    )

    // Kraken classification (using unmapped bacterial reads in FASTQ format)
    if (params.kraken) {
        kraken_classification(
            mapped_results.unmapped_fastq,
            params.kraken_db
            params.t_kraken
        )
    }

    // Quantification
    if (params.counting) {
        mapped_results.mapped_bam
            .map { bam_file -> 
                def sample_name = bam_file.simpleName.split('_mapped_')[0]
                def sample_dir = bam_file.parent
                tuple(sample_dir, sample_name)
            }
            .set { sample_info }

        counting(
            sample_info,
            params.COUNTING_SPECIES_LIST,
            params.COUNTING_ANNOTATIONS
        )
    }
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
