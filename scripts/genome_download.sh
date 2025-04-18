#!/bin/bash

### HARDCODED INPUT ####################################
OUTDIR="genomes"
FORMAT="fasta"
DATABASE="refseq"
FNA_DIR="$OUTDIR/fna_files"

### REQUIREMENTS ########################################
# Requires ncbi-genome-download (install with: pip install ncbi-genome-download)

### SPECIES LIST ########################################
# RefSeq genome 
accession_list=(
  "GCF_000159015.1" # Actinomyces coleocanis
  "GCF_002861525.1" # Actinomyces urogenitalis
  #"GCF_001543285.1" # Aerococcus viridans  
  #"GCF_000758765.1" # Anaerococcus lactolyticus  
  #"GCF_947087635.1" # Anaerococcus tetradius  
  #"GCF_016026575.1" # Atopobium vaginae  
  #"GCF_001042595.1" # Bifidobacterium dentium  
  #"GCF_000178455.1" # Brevibacterium mcbrellneri  
  #"GCF_016766875.1" # Chryseobacterium gleum  
  #"GCF_000025225.2" # Clostridiales genomosp. BVAB3  
  #"GCF_000022905.1" # Corynebacterium aurimucosum  
  #"GCF_024453835.1" # Corynebacterium genitalium  
  #"GCF_030440595.1" # Corynebacterium glucuronolyticum  
  #"GCF_028609885.1" # Corynebacterium jeikeium  
  #"GCF_000159635.1" # Corynebacterium lipophiloflavum  
  #"GCF_024453815.1" # Corynebacterium pseudogenitalium  
  #"GCF_016728105.1" # Corynebacterium striatum  
  #"GCF_947090205.1" # Dialister microaerophilus  
  #"GCF_000393015.1" # Enterococcus faecalis  
  #"GCF_000183205.1" # Eremococcus coleocola  
  #"GCF_000005845.2" # Escherichia coli  
  #"GCF_002243135.1" # Finegoldia magna  
  #"GCF_003019295.1" # Fusobacterium nucleatum  
  #"GCF_001042655.1" # Gardnerella vaginalis  
  #"GCF_034298135.1" # Lactobacillus acidophilus  
  #"GCF_009769205.1" # Lactobacillus crispatus  
  #"GCF_000056065.1" # Lactobacillus delbrueckii bulgaricus  
  #"GCF_029961225.1" # Lactobacillus fermentum  
  #"GCF_000014425.1" # Lactobacillus gasseri  
  #"GCF_011058775.1" # Lactobacillus iners  
  #"GCF_001936235.1" # Lactobacillus jensenii  
  #"GCF_014058685.1" # Lactobacillus johnsonii  
  #"GCF_025311495.1" # Lactobacillus oris  
  #"GCF_025311495.1" # Lactobacillus reuteri  
  #"GCF_006151905.1" # Lactobacillus rhamnosus  
  #"GCF_035231985.1" # Lactobacillus salivarius  
  #"GCF_035231985.1" # Lactobacillus vaginalis  
  #"GCF_000023905.1" # Leptotrichia buccalis  
  #"GCF_000177555.1" # Megasphaera genomosp. type_1  
  #"GCF_000196535.1" # Mobiluncus curtisii  
  #"GCF_000146285.1" # Mobiluncus curtisii curtisii  
  #"GCF_000185445.1" # Mobiluncus curtisii holmesii  
  #"GCF_014204735.1" # Mobiluncus mulieris  
  #"GCF_000164135.1" # Mycobacterium parascrofulaceum  
  #"GCF_000146345.1" # Peptoniphilus duerdenii  
  #"GCF_900454715.1" # Peptoniphilus lacrimalis  
  #"GCF_000212375.1" # Porphyromonas asaccharolytica  
  #"GCF_000482365.1" # Porphyromonas uenonis  
  #"GCF_000177355.1" # Prevotella amnii  
  #"GCF_000262545.1" # Prevotella bivia  
  #"GCF_025151385.1" # Prevotella buccalis  
  #"GCF_900454955.1" # Prevotella disiens  
  #"GCF_000144405.1" # Prevotella melaninogenica  
  #"GCF_000185145.2" # Prevotella oralis  
  #"GCF_000455445.1" # Prevotella timonensis  
  #"GCF_000069965.1" # Proteus mirabilis  
  #"GCF_000164635.1" # Roseomonas cervicalis  
  #"GCF_016725645.1" # Sphingobacterium spiritivorum  
  #"GCF_000013425.1" # Staphylococcus aureus  
  #"GCF_001552035.1" # Streptococcus agalactiae  
  #"GCF_013267695.1" # Streptococcus bovis  
  #"GCF_900637075.1" # Streptococcus pseudoporcinus  
  #"GCF_008153345.1" # Treponema phagedenis  
  #"GCF_902375065.1" # Veillonella atypica
)

### CREATE OUTPUT DIRECTORIES ############################
mkdir -p "$OUTDIR"
mkdir -p "$FNA_DIR"

### DOWNLOAD LOOP ########################################
for accession in "${accession_list[@]}"; do
  echo "--- Downloading: $accession ---"
  #safe_name=$(echo "$accession" | sed 's/ /_/g') # Replace spaces with underscores for safe naming

  if ! ncbi-genome-download bacteria \
      --section "$DATABASE" \
      --formats "$FORMAT" \
      --assembly-accessions "$accession" \
      --output-folder "$OUTDIR"; then
    echo "Failed to download for $accession. Skipping..."
    continue  # Skip to the next accession if download fails
  fi
    
  # Automatically unzip downloaded files
  echo "Unzipping files for $accession..."
  find "$OUTDIR" -name "*.gz" -exec gunzip -f {} \;

  # Move .fna files to the fna_files directory
  echo "Moving .fna files to $FNA_DIR..."
  find "$OUTDIR/refseq/bacteria/$accession" -name "*.fna" -exec mv {} "$FNA_DIR" \;

  echo "Done downloading and moving $accession"
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
find "$FNA_DIR" -name "*.fna" -exec cat {} + > "$OUTDIR/all_genomes.fna"
echo "All .fna files have been concatenated into $OUTDIR/all_genomes.fna"
echo # Add an empty line for better readability

### MAKE INDEX WITH MINIMAP2 ######################
echo "Creating index with minimap2..."
minimap2 -x map-pb -d $OUTDIR/bacterial_index.mmi $OUTDIR/all_genomes.fna
echo # Add an empty line for better readability
echo "Index file can be found at $OUTDIR/bacterial_index.mmi"