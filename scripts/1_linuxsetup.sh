#!/bin/bash

# Help/usage text
usage="$(basename "$0") \n
\n
This script downloads and installs Miniconda3, and uses conda to install \n
the virtual environment from GitHub.
\n
Before running this script... \n
\n
\t Please run the following command:  \n
\t \t \$ sudo dnf update \n
\t This will ensure that the software installed will be up-to-date. \n
\n
\n
Optional arguments: \n
\t      -h | --help\t         show help text and exit \n
"

# Iterating through the input arguments with a while loop
while (( "$#" )); do
	case "$1" in
		-h|--help)
			echo -e $usage
			exit 0
			;;
	esac
done

# Changing directory to the home directory
cd ~

echo Checking if the Conda environment is already installed...
sleep 2s # Slows down script to make terminal output more readable
if [ -d ~/miniconda3/envs/dRNAseq ]; then 
	echo The Conda environment already exists, exiting script.
	exit 0
else
    echo The Conda environment does not exist, continuing with script.
fi

#########################################################################################
# DOWNLOAD MINICONDA
#########################################################################################
echo Downloading Miniconda3 installation script...
sleep 2s # Slows down script to make terminal output more readable
# Check if Miniconda is installed
if command -v conda &> /dev/null; then
	echo Miniconda3 is already installed.
else
	mkdir -p ~/miniconda3 # Create Miniconda directory where installation will occur
	
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh # Install Miniconda + script
	if [ $? -ne 0 ]; then
        echo "ERROR: Failed to download Miniconda installation script."
        exit 1
    fi

	bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 # Run the script
    if [ $? -ne 0 ]; then
        echo "ERROR: Miniconda installation failed."
        exit 1
    fi

	rm ~/miniconda3/miniconda.sh # Remove the script

    echo Setting up Miniconda3...
    sleep 2s # Slows down script to make terminal output more readable
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    hash -r # Refresh the terminal after installation
    conda config --set always_yes yes --set changeps1 yes \
        --set auto_activate_base false # Answer "yes" to all prompts
    conda update -q conda # Update conda for newer versions (if available)

    #Initialize conda for the shell if it hasn't been done already
    if ! grep -q "conda initialize" ~/.bashrc; then
        conda init
    fi

    echo Displaying information about current conda installation...
    sleep 2s # Slows down script to make terminal output more readable
    conda info -a
fi

#########################################################################################
# DOWNLOAD CONDA ENVIRONMENT
#########################################################################################
echo Finding environment.yml file in the dRNASeq repository...
sleep 2s # Slows down script to make terminal output more readable

repo_dir="dRNASeq"

# Path to the extracted YAML file
env_file="$HOME/$repo_dir/environment.yml"

# Check if the YAML file exists
if [ -f "$env_file" ]; then
    echo "Found file in $env_file. Importing the Conda environment..."
    conda env create -f "$env_file"
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to import the Conda environment from $env_file."
        exit 1
    fi
else
    echo "ERROR: environment.yml file not found in the repository."
    exit 1
fi

#########################################################################################
# END SCRIPT
#########################################################################################
echo -e Script finished! \n

echo -e Please restart the terminal session for these changes to take effect. \n

echo The dRNAseq environment can be activated using the command...
echo -e	"\t \$ conda activate dRNAseq"
echo A conda virtual environment can be deactivated using the command...
echo -e	"\t \$ conda deactivate"
