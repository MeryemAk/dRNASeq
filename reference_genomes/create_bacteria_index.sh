#!/bin/bash

### HARDCODED INPUT ####################################
OUTDIR="."
FORMAT="fasta,gff"
DATABASE="refseq"
FNA_DIR="$OUTDIR/fna_files"
GFF_DIR="$OUTDIR/gff_files"

### REQUIREMENTS ########################################
# Requires:
#   - ncbi-genome-download (pip install ncbi-genome-download)
#   - mmseqs2 (conda install bioconda::mmseqs2)
#   - gffread: (conda install bioconda::gffread)
                  
### SPECIES LIST ########################################
# RefSeq genome for bacteria
bacteria=(
  "GCF_000159015.1" # Actinomyces coleocanis --> Gleimia coleocanis
  "GCF_002861525.1" # Actinomyces urogenitalis
  "GCF_001543285.1" # Aerococcus viridans
  "GCF_000758765.1" # Anaerococcus lactolyticus
  "GCF_947087635.1" # Anaerococcus tetradius
  "GCF_016026575.1" # Atopobium vaginae --> Fannyhesea vaginae
  "GCF_001042595.1" # Bifidobacterium dentium
  "GCF_000178455.1" # Brevibacterium mcbrellneri
  "GCF_016766875.1" # Chryseobacterium gleum
  "GCF_000025225.2" # Clostridiales genomosp. BVAB3 --> Mageeibacillus indolicus
  "GCF_000022905.1" # Corynebacterium aurimucosum
  "GCF_024453835.1" # Corynebacterium genitalium
  "GCF_030440595.1" # Corynebacterium glucuronolyticum
  "GCF_028609885.1" # Corynebacterium jeikeium
  "GCF_000159635.1" # Corynebacterium lipophiloflavum
  "GCF_024453815.1" # Corynebacterium pseudogenitalium
  "GCF_016728105.1" # Corynebacterium striatum
  "GCF_947090205.1" # Dialister microaerophilus --> Dialister micraerophilus
  "GCF_000393015.1" # Enterococcus faecalis
  "GCF_000183205.1" # Eremococcus coleocola
  "GCF_000005845.2" # Escherichia coli
  "GCF_002243135.1" # Finegoldia magna
  "GCF_003019295.1" # Fusobacterium nucleatum nucleatum
  "GCF_001042655.1" # Gardnerella vaginalis
  "GCF_034298135.1" # Lactobacillus acidophilus
  "GCF_009769205.1" # Lactobacillus crispatus
  "GCF_000056065.1" # Lactobacillus delbrueckii bulgaricus
  "GCF_029961225.1" # Lactobacillus fermentum
  "GCF_000014425.1" # Lactobacillus gasseri
  "GCF_011058775.1" # Lactobacillus iners
  "GCF_001936235.1" # Lactobacillus jensenii
  "GCF_014058685.1" # Lactobacillus johnsonii
  "GCF_025311495.1" # Lactobacillus oris
  "GCF_003703885.1" # Lactobacillus reuteri
  "GCF_006151905.1" # Lactobacillus rhamnosus
  "GCF_035231985.1" # Lactobacillus salivarius
  "GCF_025311515.1" # Lactobacillus vaginalis
  "GCF_000023905.1" # Leptotrichia buccalis
  "GCF_000214495.1" # Megasphaera genomosp. type 1 --> Megasphaera lornae
  "GCF_000196535.1" # Mobiluncus curtisii
  "GCF_000146285.1" # Mobiluncus curtisii curtisii
  "GCF_000185445.1" # Mobiluncus curtisii holmesii --> Mobiluncus holmesii
  "GCF_014204735.1" # Mobiluncus mulieris
  "GCF_000164135.1" # Mycobacterium parascrofulaceum
  "GCF_000146345.1" # Peptoniphilus duerdenii
  "GCF_900454715.1" # Peptoniphilus lacrimalis
  "GCF_000212375.1" # Porphyromonas asaccharolytica
  "GCF_000482365.1" # Porphyromonas uenonis
  "GCF_000177355.1" # Prevotella amnii
  "GCF_000262545.1" # Prevotella bivia
  "GCF_025151385.1" # Prevotella buccalis
  "GCF_900454955.1" # Prevotella disiens
  "GCF_000144405.1" # Prevotella melaninogenica
  "GCF_000185145.2" # Prevotella oralis
  "GCF_000455445.1" # Prevotella timonensis
  "GCF_000069965.1" # Proteus mirabilis
  "GCF_000164635.1" # Roseomonas cervicalis
  "GCF_016725645.1" # Sphingobacterium spiritivorum
  "GCF_000013425.1" # Staphylococcus aureus
  "GCF_001552035.1" # Streptococcus agalactiae
  "GCF_900102715.1" # Streptococcus bovis --> Streptococcus equines
  "GCF_900637075.1" # Streptococcus pseudoporcinus
  "GCF_008153345.1" # Treponema phagedenis
  "GCF_902375065.1" # Veillonella atypica
)
# RefSeq genome for mammals
mammals=(

)
# RefSeq genomes for fungi
fungi=(

)


### CREATE OUTPUT DIRECTORIES ############################
mkdir -p "$OUTDIR"
mkdir -p "$FNA_DIR"
mkdir -p "$GFF_DIR"

### DOWNLOAD LOOP ########################################
for accession in "${bacteria[@]}" "${mammals[@]}" "${fungi[@]}"; do
  # Check if accession already downloaded (any .fna file in $FNA_DIR contains accession)
  if ls "$FNA_DIR"/*"$accession"*.fna 1> /dev/null 2>&1; then
    echo "Accession $accession already downloaded, skipping."
    continue
  fi
  
  echo "--- Downloading: $accession ---"

  # Determine the taxonomic group based on the accession list
  if [[ " ${bacteria[@]} " =~ " $accession " ]]; then
    tax_group="bacteria"
  elif [[ " ${mammals[@]} " =~ " $accession " ]]; then
    tax_group="vertebrate_mammalian"
  elif [[ " ${fungi[@]} " =~ " $accession " ]]; then
    tax_group="fungi"
  else
    echo "Unknown taxonomic group for $accession. Exiting..."
    exit
  fi

  if ! ncbi-genome-download "$tax_group" \
      -s "$DATABASE" \
      -F "$FORMAT" \
      -A "$accession" \
      -o "$OUTDIR"; then
    echo "ncbi-genome-download $tax_group -s $DATABASE -F $FORMAT -A $accession -o $OUTDIR"
    echo "Failed to download for $accession."
    exit  # Exit script when one accession nr failes
  else
    echo "Running: ncbi-genome-download $tax_group -s $DATABASE -F $FORMAT -A $accession -o $OUTDIR"
    echo "Download completed for $accession."
  fi
    
  # Automatically unzip downloaded files
  echo "Unzipping files for $accession..."
  find "$OUTDIR" -name "*.gz" -exec gunzip -f {} \;

  # Move .fna files to the fna_files directory
  echo "Moving .fna files to $FNA_DIR..."
  find "$OUTDIR/refseq/$tax_group/$accession" -name "*.fna" -exec mv {} "$FNA_DIR" \;

  # Move .gff files to the gff_files directory
  echo "Moving .gff files to $GFF_DIR..."
  find "$OUTDIR/refseq/$tax_group/$accession" -name "*.gff" -exec mv {} "$GFF_DIR" \;

  echo "Done with $accession"
  echo # Add an empty line for better readability
done

### REMOVE UNNECESSARY FILES ############################
echo "Removing unnecessary files..."
chmod -R u+w "$OUTDIR/refseq"
rm -rf "$OUTDIR/refseq"
echo "Unnecessary files removed."
echo # Add an empty line for better readability

### CONCATENATE ALL .fna FILES INTO ONE ###################
echo "Concatenating all .fna files into one file..."
if [ ! -f "$OUTDIR/bacteria_seq.fna" ]; then
  find "$FNA_DIR" -name "*.fna" -exec cat {} + > "$OUTDIR/bacteria_seq.fna"
  echo "All .fna files have been concatenated into $OUTDIR/bacteria_seq.fna"
else
  echo "Concatenated .fna file already exists. Skipping this step."
fi
echo # Add an empty line for better readability

### CLUSTER SEQUENCES WITH MMseqs2 ###################
echo "Clustering sequences with MMseqs2..."
if [ ! -f "$OUTDIR/bacteria_seq_mmseqs_rep_seq.fasta" ]; then
  mmseqs easy-linclust "$OUTDIR/bacteria_seq.fna" "$OUTDIR/bacteria_seq_mmseqs" tmp --min-seq-id 0.97 --cov-mode 1 -c 0.97 --split-memory-limit 4000 --threads 2
  echo "Clustering complete. Output: $OUTDIR/bacteria_seq_mmseqs_rep_seq.fasta"
else
  echo "MMseqs output already exists. Skipping clustering step."
fi
echo # Add an empty line for better readability

### MAKE INDEX WITH MINIMAP2 ######################
echo "Creating index with minimap2..."
minimap2 -x map-ont -d $OUTDIR/bacteria_index.mmi $OUTDIR/bacteria_seq_mmseqs_rep_seq.fasta
echo # Add an empty line for better readability
echo "Index file can be found at $OUTDIR/bacteria_index.mmi"

### CONVERT GFF FILES TO GTF FORMAT ##################
echo "Converting .gff files to .gtf format..."
if [ ! -f "$OUTDIR/bacteria_annotation.gtf" ]; then
  mkdir -p "$OUTDIR/gtf_files"
  for gff_file in "$GFF_DIR"/*.gff; do
    gffread "$gff_file" -T -o "$OUTDIR/gtf_files/$(basename "${gff_file%.gff}.gtf")"
  done
  echo "Converted GFF files to GTF format in $OUTDIR/gtf_files"
else
  echo "GTF files already exist. Skipping conversion step."
fi

### CONCATENATE ALL .gtf FILES INTO ONE ###################
echo "Concatenating .gtf files into one..."
if [ ! -f "$OUTDIR/bacteria_annotation.gtf" ]; then
  find "$GFF_DIR" -name "*.gtf" -exec cat {} + > "$OUTDIR/bacteria_annotation.gtf"
  echo "Annotation file: $OUTDIR/bacteria_annotation.gtf"
else
  echo "Concatenated .gtf file already exists. Skipping this step."
fi
echo # Add an empty line for better readability
