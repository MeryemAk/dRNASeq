#!/usr/bin/env python3
#########################################################################################
# Import libraries
import gzip
import os
#########################################################################################
# HARDCODE PARAMETERS
baseDir = os.getcwd()                          # Get the current directory
inFolder = os.path.join(baseDir, "sequences")  # Define the sequences folder path
outFolder = "trimming"
#########################################################################################
# CHECK FOLDERS
# Check if inFolder exists
if not os.path.exists(inFolder):
    print(f"Error: The folder '{inFolder}' does not exist.")
    exit(1)  # Exit the script with an error code

os.makedirs(outFolder, exist_ok=True) # Create the output folder if it doesn't exist
#########################################################################################
# FILENAMES
fileNames = []
for file in os.listdir(inFolder): # get files with fastq.gz extension
    if file.endswith(".fastq.gz"):
        fileNames.append(os.path.join(inFolder, file)) # Append the file name to the list
        print("Found files:")
        print("--- {} ---".format(file))
    else:
        print(f"Error: No .fastq.gz files found in '{inFolder}'.")
        exit(1)
#########################################################################################
# RUN FILTLONG --> LONG READS
# Change parameters if necessary
min_length = 1000           # Minimum read length                   # --min_length 1000 ← Discard any read shorter than 1 kbp.
keep_percent = 90           # Percentage of best reads to keep      #--keep_percent 90 ← Throw out the worst 10% of reads. (measured by bp not read count)
target_bases = 500000000    # Target number of bases                # --target_bases 500000000 ← Remove the worst reads until only 500 Mbp remain, useful for very large read sets.

# Execute command for each file
for file in fileNames:
    # Generate output filepath inside outFolder
    output_file = os.path.join(outFolder, os.path.basename(file))

    # Compose command
    filt_cmd = (
        f"filtlong --min_length {min_length} "
        f"--keep_percent {keep_percent} "
        f"--target_bases {target_bases} "
        f"{file} | gzip > {output_file}"
    )

    print("\nexectuing: {}".format(filt_cmd))
    exit_code = os.system(filt_cmd) # Execute

    # Check if the command was successful
    if exit_code != 0:
        print(f"Error: Filtlong failed with exit code {exit_code}.")

print("Filtlong analysis completed.\n")
