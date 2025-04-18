#!/usr/bin/env python3
#########################################################################################
# Import libraries
import os      # file and directory manipulation
import time    # time manipulation
#########################################################################################
baseDir = os.getcwd()                          # Get the current directory
inFolder = os.path.join(baseDir, "sequences")  # Define the sequences folder path
threads = 4                                    # Number of threads to use 
#########################################################################################
# INPUT FOLDER
#########################################################################################
# Check if inFolder exists
if not os.path.exists(inFolder):
    print(f"Error: The folder '{inFolder}' does not exist.")
    exit(1)  # Exit the script with an error code
#########################################################################################
# OUTPUT FOLDER
#########################################################################################
outFolder = "QualityControl"
os.makedirs(outFolder, exist_ok=True) # Create the output folder if it doesn't exist
#########################################################################################
# FILENAMES
#########################################################################################
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
# RUN FASTQC --> SHORT READS
#########################################################################################
# Compose command
for file in fileNames:
    fastqc_cmd = "fastqc --extract -t {} -o {} {}".format(threads, outFolder, file)
    print("\nexectuing: {}".format(fastqc_cmd))
    exit_code = os.system(fastqc_cmd) # Execute

    # Check if the command was successful
    if exit_code != 0:
        print(f"Error: FastQC failed with exit code {exit_code}.")
    else:
        print(f"FastQC analysis completed successfully for {file}.\n")

print("FastQC analysis completed.\n")

#########################################################################################
# RUN LONGQC --> LONG READS
#########################################################################################
# Construct LongQC command
#for file in fileNames:
#    longqc_cmd = f"longQC analysis standard --out_dir {outFolder} --threads {threads} {file}"
#    print(f"\nExecuting: {longqc_cmd}")
#    exit_code = os.system(longqc_cmd)  # Execute LongQC
#
#    # Check if the command was successful
#    if exit_code != 0:
#        print(f"Error: LongQC failed for {file} with exit code {exit_code}.")
#    else:
#        print(f"LongQC analysis completed successfully for {file}.\n")

#print("LongQC analysis for all files completed.")

#########################################################################################
# RUN MULTIQC
#########################################################################################
print("Creating MultiQC report...\n")

os.chdir(outFolder)  # Change to the output folder for MultiQC
print("Moving into: {}".format(outFolder))

multiqc_cmd = "multiqc ."
os.system(multiqc_cmd)  # Run MultiQC in the output folder

print("MultiQC report created.")

print("Deleting unnecessary files...")
time.sleep(2)
delete_cmd="rm *.zip"
os.system(delete_cmd)

print("The script is done!")

