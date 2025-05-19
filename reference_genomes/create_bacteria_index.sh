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

### SPECIES LIST ########################################
# RefSeq genome for bacteria
bacteria=(
  "GCF_000005845.2" # Escherichia coli
  "GCF_000013425.1" # Staphylococcus aureus
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
find "$FNA_DIR" -name "*.fna" -exec cat {} + > "$OUTDIR/bacteria_seq.fna"
echo "All .fna files have been concatenated into $OUTDIR/bacteria_seq.fna"
echo # Add an empty line for better readability

### MAKE INDEX WITH MINIMAP2 ######################
echo "Creating index with minimap2..."
minimap2 -x map-pb -d $OUTDIR/bacteria_index.mmi $OUTDIR/bacteria_seq.fna
echo # Add an empty line for better readability
echo "Index file can be found at $OUTDIR/bacteria_index.mmi"

### CONCATENATE ALL .gff FILES INTO ONE ###################
echo "Concatenating .gff files into one..."
find "$GFF_DIR" -name "*.gff" -exec cat {} + > "$OUTDIR/bacteria_annotation.gff"
echo "Annotation file: $OUTDIR/bacteria_annotation.gff"
