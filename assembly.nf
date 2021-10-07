#!/usr/bin/env nextflow

/* Help message*/
if (params.help) {
    print helpMessage()
    exit 0
}

else {
    print parameterShow()
}

/* Parameters validation */
error = false
if (!params.in_dir){
    println "${c_red}\nParameter --in_dir is missing${c_reset}"
    error = true
} 

else if (params.in_dir instanceof Boolean) {
    println "${c_red}\nParameter --in_dir is empty, no input directory specified${c_reset}"
    error = true
} 

else if (!file(params.in_dir).exists()){
    println "${c_red}\nParameter --in_dir is not an existing directory${c_reset}"
    error = true
}

if (!params.out_dir){
    println "${c_red}\nParameter --out_dir is missing${c_reset}"
    error = true
} 

else if (params.out_dir instanceof Boolean) {
    println "${c_red}\nParameter --out_dir is empty, no output directory specified${c_reset}"
    error = true
} 

else if (!file(params.out_dir).isDirectory()){
    println "${c_red}\nParameter --out_dir is not an existing directory${c_reset}"
    error = true
}

if (params.basecalling && (file("${params.in_dir}/fast5/").isDirectory())) {
    count = 0
    file("${params.in_dir}/fast5").eachFileRecurse {
        if (it.name.endsWith(".fast5")) {
            count += 1
        }
    }
    if (count == 0) {
        println "${c_red}\nNo fast5 files found in the fast5 directory${c_reset}"
        error=true
    }
}

if (!params.basecalling && !params.merged && (!file("${params.in_dir}/fastq/").isDirectory())) {
    println "${c_red}\nParameter --in_dir needs to contain a fastq directory with fastq files${c_reset}"
    println "${c_red}Make a fastq directory inside the in_dir${c_reset}"
    error=true
}

if (!params.basecalling && !params.merged && (file("${params.in_dir}/fastq/").isDirectory())) {
    count = 0
    file("${params.in_dir}/fastq").eachFileRecurse {
        if (it.name.endsWith(".fastq")) {
            count += 1
        }
    }
    if (count == 0) {
        println "${c_red}\nNo fastq files found in the fastq directory${c_reset}"
        error=true
    }
}

if (!params.basecalling && !params.merged && (file("${params.in_dir}/fastq/merged").isDirectory())) {
    println "${c_red}\nParameter --in_dir contains a fastq/merged directory, are fastq's already merged? ${c_reset}"
    println "${c_red}Set parameter --merged to true ${c_reset}"
    error=true
}

if (!params.basecalling && params.merged && (!file("${params.in_dir}/fastq/merged").isDirectory())) {
    println "${c_red}\nParameter --in_dir needs to contain a fastq/merged directory when parameter --merged true${c_reset}"
    error=true
}

else if (!params.basecalling && params.merged && (file("${params.in_dir}/fastq/merged").isDirectory())) {
    count = 0
    file("${params.in_dir}/fastq/merged").eachFile {
        if (it.name.endsWith(".fastq")) {
            count += 1
        }
    }
    if (count == 0) {
        println "${c_red}\nNo fastq files found in the fastq/merged directory${c_reset}"
        error=true
    }
}

if (params.asm_coverage && !params.gsize){
    println "${c_red}\nIf params.asm_coverage is specified, params.gsize is mandatory${c_reset}"
    error=true
}

if (params.polishing && !params.mapping){
    println "${c_red}\nIf parameter --polishing is provided, parameter --mapping is mandatory${c_reset}"
    error=true
}

if (params.miniasm_assembly && params.flye_assembly){
    println "${c_red}\nTwo assembly methods are specified, please choose only one${c_reset}"
    error=true
}

if (error) {
    exit 0
}


/* processes */
process guppy {
    
    label "guppy"
    
    publishDir "${params.out_dir}/01.basecalled/", mode: "copy"
    
    input:
    path x from Channel.fromPath(params.in_dir+"/fast5")
    
    output:
    path "merged/*.fastq" into fastq_w_guppy
    path "*" 

    when:
    params.basecalling

    script:
    barcode_kits = params.barcode_kits ? "--barcode_kits $params.barcode_kits" : ""
    skip_qc = params.skip_qc ? "" : "--qscore_filtering"

    """
    IFS=','
    read -ra BCS <<<"${params.barcodes}"
    # 1 barcode with leading 0: octal for bash
    if [[ \${#BCS[@]} -eq 1 ]] && [[ \${#BCS[0]} -eq 1 ]]; then
        BCS[0]="0"\${BCS[0]}
    fi
    
    for d in ./fast5/*
    do
        # guppy without barcodes
        if [[ -f \$d && \$d == *.fast5 ]]; then
            guppy_basecaller \
            --input_path ${params.in_dir}/fast5 \
            --save_path ./ \
            -c ${params.bc_config} \
            --num_callers ${params.num_callers} \
            $barcode_kits \
            $skip_qc 
            break
        # guppy with barcodes
        elif [[ \$d =~ barcode* ]]; then
            # all the barcodes
            if [[ \${#BCS[@]} -eq 0 ]]; then
                guppy_basecaller \
                --input_path \$d \
                --save_path ./ \
                -c ${params.bc_config} \
                --num_callers ${params.num_callers} \
                $barcode_kits \
                $skip_qc 
            # guppy with barcodes from params.barcodes
            else
                bc=\${d##*barcode}
                if [[ " \${BCS[@]} " =~ " \$bc " ]]; then
                    guppy_basecaller \
                    --input_path \$d \
                    --save_path ./ \
                    -c ${params.bc_config} \
                    --num_callers ${params.num_callers} \
                    $barcode_kits \
                    $skip_qc 

                fi
            fi
        fi
    done

    # merge the fastq files  
    mkdir merged

    if [[ ${params.skip_qc} == false ]]; then
        directory="./pass/*"
    else
        directory="./*"
    fi

    for dir in \$directory
    do
        if [[ -f \$dir && \$dir == *.fastq ]]; then
            cat ./*.fastq > ./merged/sample.fastq
            break
            
        elif [[ \$dir =~ barcode* ]]; then
            cat \$dir"/"*.fastq > ./merged/\$(basename \$dir)".fastq"
        fi
    done
    """
}

process merge {

    label "merge"

    publishDir "${params.out_dir}/01.basecalled/", mode: "copy"

    input: 
    path fq from Channel.fromPath(params.in_dir+"/fastq").collect()

    when:
    !params.merged && !params.basecalling

    output:
    path "merged/*.fastq" into fastq_wo_guppy

    """
    mkdir merged
    IFS=','
    read -ra BCS <<<"${params.barcodes}"
     # 1 barcode with leading 0: octal for bash
    if [[ \${#BCS[@]} -eq 1 ]] && [[ \${#BCS[0]} -eq 1 ]]; then
        BCS[0]="0"\${BCS[0]}
    fi

    for d in ./fastq/*
    do
        # merge without barcodes
        if [[ -f \$d ]]; then
            cat ./fastq/*.fastq > ./merged/sample.fastq
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

    label "no_merge"

    publishDir "${params.out_dir}/01.basecalled/merged", mode: "copy"
    
    input: 
    path fq from Channel.fromPath(params.in_dir+"/fastq/merged/*")

    when:
    params.merged && !params.basecalling

    output:
    path "**.fastq" optional true into fastq_merged

    """
    IFS=','
    read -ra BCS <<<"${params.barcodes}"
     # 1 barcode with leading 0: octal for bash
    if [[ \${#BCS[@]} -eq 1 ]] && [[ \${#BCS[0]} -eq 1 ]]; then
        BCS[0]="0"\${BCS[0]}
    fi

    mkdir out
    dir="./*"
    # make directory to copy the merged fastq 
    # when filenames are same --> no problem copying them to other directory
    
    for fastq in \$dir
    do
        if [[ \${#BCS[@]} -ge 1 ]]; then
            for i in "\${BCS[@]}"
            do
                if [[ "\${fastq}" =~ "\${i}" ]]
                then 
                    cat \${fastq} > ./out/\${fastq} 
                fi
            done
        # als er geen barcodes meegegeven worden
        elif [[ \${#BCS[@]} -eq 0 ]] && [[ ! -d \$fastq ]]; then
            cat \${fastq} > ./out/\${fastq} 
        fi
    done
    """
}

fastqs = fastq_w_guppy.mix(fastq_wo_guppy, fastq_merged)
fastqs.into{fastq_nanoplot;fastq_nanocomp;fastq_flye;fastq_miniasm}

process nanoplot {

    label "nanoplot"

    /* Ignore error for limited data*/
    errorStrategy "ignore"

    publishDir "${params.out_dir}/02.qc/nanoplot/${fq.simpleName}", mode:"copy", pattern: "*.png"
    publishDir "${params.out_dir}/02.qc/nanoplot/${fq.simpleName}", mode:"copy", pattern: "*.txt"
    publishDir "${params.out_dir}/02.qc/nanoplot/${fq.simpleName}", mode:"copy", pattern: "*.html"

    input:
    // .flatten() analyses all the barcodes seperate
    path fq from fastq_nanoplot.flatten()

    output:
    path "*" optional true
    tuple file(fq), file("*.txt"), file("HistogramReadlength.png"), file("LengthvsQualityScatterPlot_kde.png") optional true into nanoplot_report

    when:
    params.nanoplot

    script:
    """
    NanoPlot -t ${params.t_qc} -o . --N50 --fastq_rich $fq
    """
}

process nanocomp {

    label "nanocomp"

    publishDir "${params.out_dir}/02.qc/nanocomp/", mode:"copy"

    input:
    // .collect() analyses all the barcodes together
    path fq from fastq_nanocomp.collect()

    output:
    path "*"
    
    when:
    params.nanocomp 

    script:
    """
    NanoComp -t ${params.t_qc} -o . -f png --plot violin --fastq $fq
    """
}

process flye_assembly {
        
    label "flye"

    /* when data is limited, sometimes Flye produces an assembly empty assembly, which will cause an error.
    errorStragy 'ignore' will not terminate the script if this occurs */
    errorStrategy "ignore"
        
    publishDir "${params.out_dir}/03.assembly/", mode: "copy", pattern: "**.fasta"
    publishDir "${params.out_dir}/03.assembly/", mode: "copy", pattern: "**.gfa"
    publishDir "${params.out_dir}/03.assembly/", mode: "copy", pattern: "**.gv"
    publishDir "${params.out_dir}/03.assembly/", mode: "copy", pattern: "**.txt"
    publishDir "${params.out_dir}/03.assembly/", mode: "copy", pattern: "**.log"

    input: 
    path fq from fastq_flye.flatten()

    output:
    path "*/*assembly*" optional true
    path "**.log" optional true
    tuple file(fq), file("*/*assembly.fasta") optional true into (flye_prokka, flye_mapping)

    when:
    params.flye_assembly

    script:

    gsize = params.gsize ? "--genome-size $params.gsize" : ""
    plasmids = params.plasmids ? "--plasmids" : ""
    coverage = params.asm_coverage ? "--asm-coverage $params.asm_coverage" : ""
    meta = params.meta ? "--meta" : ""

    """
    flye --nano-raw $fq $gsize $coverage $plasmids $meta --out-dir ./${fq.simpleName} -t ${params.t_assembly} -i 5

    # add bc on assembly.fasta
	for file in ./${fq.simpleName}/assembly*; do
		basename=\${file##*/}
		basename_no_ext=\${basename%.*}
		extension=\${file##*.}
		mv "\$file" "./${fq.simpleName}/${fq.simpleName}_\$basename_no_ext.\$extension"
	done
    """
} 

process miniasm_assembly {

    label "miniasm"

    publishDir "${params.out_dir}/03.assembly/", mode: "copy", pattern: "**.gfa"
    publishDir "${params.out_dir}/03.assembly/", mode: "copy", pattern: "**.fasta"

    input:
    path fq from fastq_miniasm.flatten()

    output:
    path "**" into miniasm_assembly
    tuple file(fq), file("**/*.fasta") optional true into (miniasm_prokka, miniasm_mapping)

    when:
    params.miniasm_assembly

    shell:
    '''
    mkdir ./!{fq.simpleName}
    minimap2 -x ava-ont -t!{params.t_assembly} !{fq} !{fq} | gzip -1 > ./!{fq.simpleName}/!{fq.simpleName}_reads.paf.gz
    miniasm -f !{fq} ./!{fq.simpleName}/!{fq.simpleName}_reads.paf.gz > ./!{fq.simpleName}/!{fq.simpleName}.gfa
    awk '/^S/{print ">"$2"\\n"$3}' ./!{fq.simpleName}/!{fq.simpleName}.gfa > ./!{fq.simpleName}/!{fq.simpleName}.fasta
    '''
}

if (params.mapping && params.flye_assembly) {
    mapping_tuple = flye_mapping
}
else if (params.mapping && params.miniasm_assembly) {
    mapping_tuple = miniasm_mapping
}
else {
    mapping_tuple = tuple()
}

process mapping {

    label "mapping"

    publishDir "${params.out_dir}/04.mapping/", mode: "copy", pattern: "*.sam"

    input:
    tuple file(fq), file(assembly) from mapping_tuple

    output:
    path "*" into sam_files
    tuple file(fq), file(assembly), file("*.sam") into racon

    when:
    params.mapping || params.polishing
    
    """
    minimap2 -a $assembly $fq -t ${params.t_mapping} > ${fq.simpleName}.sam
    """
}

process sam_to_bam {

    label "samtools"

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

    label "samtools"

    publishDir "${params.out_dir}/04.mapping", mode: "copy"

    input:
    path x from bam_files

    output:
    path "*" into sorted_bam_files

    when:
    params.mapping
    
    """
    samtools sort $x -o ${x.simpleName}_sorted
    """
}

process index_bam {
    
    label "samtools"

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

    label "racon"
    
    errorStrategy "ignore"
   
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
    racon $fq $alignment $assembly -t ${params.t_polishing} -m 8 -x -6 -g -8 -w 500 > ${fq.simpleName}_racon.fasta 
    """
}

process medaka {

    label "medaka"

    publishDir "${params.out_dir}/05.polishing/medaka", mode: "copy", pattern: "**.fasta"

    input:
    tuple file(fq), file(assembly) from medaka

    output:
    path "**.fasta" into medaka_result
    tuple file(fq), file ("**consensus.fasta") into medaka_prokka

    when:
    params.polishing

    """
    medaka_consensus -i $fq -d $assembly -o ./${fq.simpleName} -t ${params.t_polishing} -m ${params.model} 
    mv ./${fq.simpleName}/consensus.fasta ${fq.simpleName}_consensus.fasta
    """

}

if (params.polishing) {
    prokka_tuple = medaka_prokka
}
else if (!params.polishing && params.flye_assembly) {
    prokka_tuple = flye_prokka
}
else if (!params.polishing && params.miniasm_assembly) {
    prokka_tuple = miniasm_prokka
}

process prokka_annotation {
  
    label "prokka"
    
    publishDir "${params.out_dir}/06.annotation/", mode: "copy"

    input:
    tuple file(fq), file(assembly) from prokka_tuple

    output:
    path "*" into prokka_result

    when:
    params.annotation

    """
    prokka --outdir ./${fq.simpleName} --prefix ${fq.simpleName} $assembly --cpus ${params.t_annotation}
    """
}

process report {

    label "report"

    publishDir "${params.out_dir}/07.report", mode: "copy"

    input:
    tuple file(fq), file(summary), file(histogram), file(quality) from nanoplot_report

    output:
    path "*" into report

    when:
    params.nanoplot

    """
    #!/usr/bin/env python

    # import modules
    from fpdf import FPDF

    class PDF(FPDF):
        def head_title(self, head_title):
            # set style
            self.set_font('Arial', 'B', 18)
            self.ln(6)
            self.cell(w=0,txt = head_title)
            self.ln(6)

        def chapter_title(self, title):
            # set style
            self.set_font('Arial', 'B', 14)
            self.ln(10)
            self.cell(w=0,txt = title)
            # Line break
            self.ln(6)

        def table(self, var, value):
            # set style
            self.set_font('Arial', '', 12)
            # cell(w, h = 0, txt = '', border = 0, ln = 0, align = '', fill = False, link = '')
            self.cell(50,10,txt = var, border = 1, ln= 0)
            self.cell(50,10,txt = value, border = 1, ln= 0)
            self.ln(10)

    with open("$summary", "r") as f:
        lines = f.readlines()
        mean_read_length = lines[2].split()[3]
        mean_read_quality = lines[3].split()[3]
        number_of_reads  = lines[6].split()[3]
        read_length_N50  = lines[7].split()[3]
        total_bases = lines[8].split()[2]
        longest_reads_1 = lines[22].split()[1] + lines[22].split()[2] 
        longest_reads_2 = lines[23].split()[1] + lines[23].split()[2] 
        longest_reads_3 = lines[24].split()[1] + lines[24].split()[2] 
        longest_reads_4 = lines[25].split()[1] + lines[25].split()[2] 
        longest_reads_5 = lines[26].split()[1] + lines[26].split()[2] 
    
    pdf = PDF()
    pdf.add_page()
    pdf.head_title("${fq.simpleName}")


    pdf.chapter_title("Tools used")
    if ("$params.basecalling" == "true"):
        pdf.table("Basecalling", "Guppy")
    elif ("$params.basecalling" == "false"):
        pdf.table("Basecalling", "None")

    if ("$params.nanoplot" == "true") and ("$params.nanocomp" == "true"):
        pdf.table("QC", "NanoPlot + NanoComp")
    elif ("$params.nanoplot" == "true") and ("$params.nanocomp" == "false"):
        pdf.table("QC", "NanoPlot")
    
    if ("$params.miniasm_assembly" == "true"):
        pdf.table("Assembly", "Miniasm")
    elif ("$params.flye_assembly" == "true"):
        pdf.table("Assembly", "Flye")
    elif ("$params.flye_assembly" == "false") and ("$params.flye_assembly" == "false"):
        pdf.table("Assembly", "None")
    
    if ("$params.mapping" == "true") and (("$params.flye_assembly" == "true") or ("$params.miniasm_assembly" == "true")) :
        pdf.table("Mapping", "Minimap2 + Samtools")
    else:
        pdf.table("Mapping", "None")
    
    if ("$params.polishing" == "true" and (("$params.flye_assembly" == "true") or ("$params.miniasm_assembly" == "true")) ):
        pdf.table("Polishing", "Medaka + Racon")
    else:
        pdf.table("Polishing", "None")

    if ("$params.annotation" == "true") and (("$params.flye_assembly" == "true") or ("$params.miniasm_assembly" == "true")) :
        pdf.table("Annotation", "Prokka")
    else:
        pdf.table("Annotation", "None")


    pdf.chapter_title("General statistics")
    pdf.table("Mean read length", mean_read_length)
    pdf.table("Mean read quality", mean_read_quality)
    pdf.table("Number of reads", number_of_reads)
    pdf.table("Read length N50", read_length_N50)
    pdf.table("Total bases", total_bases)
    pdf.chapter_title("Top 5 longest reads and their mean basecall quality score")
    pdf.table("1", longest_reads_1)
    pdf.table("2", longest_reads_2)
    pdf.table("3", longest_reads_3)
    pdf.table("4", longest_reads_4)
    pdf.table("5", longest_reads_5)

    pdf.add_page()
    pdf.chapter_title("Figures")
    pdf.image("$histogram", w = 175)
    pdf.image("$quality", w = 175)


    pdf.output('${fq.simpleName}.pdf', 'F')
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
    • in_dir:                       | ${c_yellow}Input directory (directory with fast5 or fastq folder)${c_reset}
    • out_dir:                      | ${c_yellow}Output directory with results${c_reset}

                                                       | ${c_blue}Optional parameter${c_reset}
    • barcodes:                     | ${c_yellow}Comma separated list of barcode numbers that the user wants to analyse if barcodes are present${c_reset}
                                    | ${c_yellow}All barcodes are automatically analysed if barcodes are present but no --barcodes option is provided${c_reset}
                                    | ${c_yellow}Numbers should include the leading 0s. E.g. 03,08,11${c_reset}

    ${c_cyan}Basecalling:${c_reset}                    | ${c_blue}Only when providing --basecalling option${c_reset}
    • basecalling:                  | ${c_yellow}If provided, this will basecall fast5 files from in_dir ${c_reset}
   	default: false
    • barcode_kits:                 | ${c_yellow}If provided, this will demultiplex the samples using the provided kit${c_reset}
    • bc_config:                    | ${c_yellow}If provided, will use a different config file then default${c_reset}
        default: dna_r9.4.1_450bps_fast.cfg
    • skip_qc:                      | ${c_yellow}If provided, this will not split basecalled reads into pass or fail${c_reset}
   	default: false 
    • num_callers (default:1):      | ${c_yellow}Number of callers to use for basecalling (1=4 threads!)${c_reset}
	default: 1

    ${c_cyan}Merging:${c_reset}  
    • merged:                       | ${c_yellow}If provided: will expect a /fastq folder with fastq files for each barcode / one sample${c_reset}
    	default: false              | ${c_yellow}If not provided: will expect a /fastq/merged/ folder with the merged fastq files for each barcode${c_reset}
                                    | ${c_yellow}/ a sample${c_reset}
   
    ${c_cyan}QC:${c_reset}                             | ${c_blue}Quality Control steps${c_reset}
    • nanocomp:                     | ${c_yellow}If provided, will perform nanocomp analysis ${c_reset}
	default: true  
    • nanoplot:                     | ${c_yellow}If provided, will perform nanoplot analysis ${c_reset}
	default: true
    • t_qc      	 	            | ${c_yellow}Number of threads used for QC${c_reset}

    ${c_cyan}assembly:${c_reset}                       | ${c_blue}Asssembly option${c_reset}
    • flye_assembly:                | ${c_yellow}If provided, this will assemble the genomes using Flye ${c_reset}
	default:true
    • gsize:                        | ${c_yellow}Expected genome size (not mandatory) ${c_reset}
    • meta:                         | ${c_yellow}Metagenome / Uneven coverage${c_reset}
    • plasmids:                     | ${c_yellow}Rescue short unassembled plasmids${c_reset}
    • asm_coverage:                 | ${c_yellow}Reduced coverage for initial disjointig assembl${c_reset}
    • miniasm_assembly:             | ${c_yellow}If provided, this will assemble the genomes using Miniasm ${c_reset}
	default:false
    • t_assembly	 	            | ${c_yellow}Number of threads per barcode (use max: 32/nbarcodes) ${c_reset}
	default: 4   


    ${c_cyan}mapping:${c_reset}                        | ${c_blue}If provided, will map, sort and index sequences ${c_reset}
	default: true
    • t_mapping                     | ${c_yellow}Number of threads used for mapping ${c_reset}

    ${c_cyan}polishing:${c_reset}                      | ${c_blue}If provided, will polish sequences (requires mapping) ${c_reset}
	default: true   
    • t_polishing                   | ${c_yellow}Number of threads used for polishing ${c_reset}
   	default: 4
    • model:                        | ${c_yellow}Model used for Medaka polishing: {pore}_{device}_{caller variant}_{caller version} ${c_reset}
        default: r941_min_fast_g303
	    hac: r941_min_high_g360

    ${c_cyan}annotation:${c_reset}                     | ${c_blue}If provided, will anotate sequences ${c_reset}
   	default: true
    • t_annotation                  | ${c_yellow}Number of threads used for annotation ${c_reset} 
        default: 4 
 
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
    • barcodes:                     | ${c_yellow}${params.barcodes} ${c_reset}

    ${c_cyan}Basecalling:${c_reset} 
    • basecalling:                  |"""

    if (params.basecalling) {
        text = text + """ ${c_yellow}${params.basecalling}${c_reset}
    • barcode_kits:                 | ${c_yellow}${params.barcode_kits}${c_reset}
    • bc_config:                    | ${c_yellow}${params.bc_config}${c_reset}
    • skip_qc:                      | ${c_yellow}${params.skip_qc}${c_reset}
    • num_callers:                  | ${c_yellow}${params.num_callers}${c_reset}

    """
    }
    else {
        text = text+""" ${c_yellow}${params.basecalling}${c_reset}

    """
    }

    text = text + """${c_cyan}Merge:${c_reset}
    • merged:                       | ${c_yellow}${params.merged}${c_reset}

    """

    text = text + """${c_cyan}QC:${c_reset}
    • nanocomp:                     | ${c_yellow}${params.nanocomp}${c_reset}
    • nanoplot:                     | ${c_yellow}${params.nanoplot}${c_reset}
    • t_assembly:                   | ${c_yellow}${params.t_qc} ${c_reset}

    ${c_cyan}Assembly:${c_reset}
    • flye_assembly:                |"""
    if (params.flye_assembly) {
        text = text + """ ${c_yellow}${params.flye_assembly}${c_reset}
    • miniasm_assembly:             | ${c_yellow}${params.miniasm_assembly}${c_reset}
    • gsize:                        | ${c_yellow}${params.gsize}${c_reset}
    • meta:                         | ${c_yellow}${params.meta}${c_reset}
    • plasmids:                     | ${c_yellow}${params.plasmids}${c_reset}
    • asm_coverage:                 | ${c_yellow}${params.asm_coverage}${c_reset}
    • t_assembly:                   | ${c_yellow}${params.t_assembly} ${c_reset}
    """
    }

    else if (params.miniasm_assembly) {
    text = text + """ ${c_yellow}${params.flye_assembly}${c_reset}
    • miniasm_assembly:             | ${c_yellow}${params.miniasm_assembly}${c_reset}
    • t_assembly:                   | ${c_yellow}${params.t_assembly} ${c_reset}  
    """
    } 
    else {
        text = text+""" ${c_yellow}${params.miniasm_assembly}${c_reset}

    """
    }
    
    text = text + """ ${c_cyan}Mapping:${c_reset}
    • mapping:                      |"""

    if (params.mapping) {
        text = text + """ ${c_yellow}${params.mapping}${c_reset}
    • t_mapping :                   | ${c_yellow}${params.t_mapping} ${c_reset}

    """
    }
    else {
        text = text+""" ${c_yellow}${params.mapping}${c_reset}

    """
    }

    text = text + """${c_cyan}Polishing:${c_reset}
    • polishing:                    |"""

    if (params.polishing) {
        text = text + """ ${c_yellow}${params.polishing}${c_reset}
    • model:                        | ${c_yellow}${params.model} ${c_reset}
    • t_polishing:                  | ${c_yellow}${params.t_polishing} ${c_reset}

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
    • t_annotation:                 | ${c_yellow}${params.t_annotation}${c_reset} 

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

