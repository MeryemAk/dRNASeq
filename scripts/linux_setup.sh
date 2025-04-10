#! /bin/bash/

# Help/usage text
usage="$(basename "$0") \n
\n
This script downloads and installs Miniconda3, and uses conda to install \n
the virtual environment from GitHub.
\n
Before running this script... \n
\n
\t 1. Please run the following command:  \n
\t \t \$ sudo dnf update \n
\t This will ensure that the software installed will be up-to-date. \n
\n
\t 2. Ensure the GitHub environment directory is in the home directory. \n
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

echo Checking if the GitHub environment is already installed...
sleep 2s # Slows down script to make terminal output more readable
if [ -d ~/miniconda/envs/githubenv ]; then 
	echo The GitHub environment already exists, exiting script.
	exit 0
fi

echo Checking if GitHub environment  is in the home directory...
sleep 2s # Slows down script to make terminal output more readable
# If githubenv/ is not in the home directory...
if [ ! -d ~/githubenv/ ];
then
	echo ERROR: GitHub environment is not in the home directory
	echo The home directory is $HOME
	echo Please move the GitHub environment directory to the home directory,
	echo or create a copy of the GitHub environment in $HOME
	exit 1
fi

echo Downloading Miniconda3 installation script...
sleep 2s # Slows down script to make terminal output more readable
# Check if Miniconda is installed
if command -v conda &> /dev/null; then
	echo Miniconda3 is already installed, exiting script.
	exit 0
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
fi

echo Setting up Miniconda3...
sleep 2s # Slows down script to make terminal output more readable
source "$HOME/miniconda/etc/profile.d/conda.sh"
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

echo Creating the GitHub virtual environment using conda...
sleep 2s # Slows down script to make terminal output more readable
if [ -f ~/path/to/GitHub/environment ]; then
    conda create --name environment_name --file ~/path/to/GitHub/environment
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to create the Conda environment."
        exit 1
    fi
else
    echo "ERROR: Environment file not found at ~/path/to/GitHub/environment"
    exit 1
fi

echo Removing unused packages and caches using conda...
sleep 2s # Slows down script to make terminal output more readable
conda clean --all --yes

echo -e Script finished! \n

echo -e Please restart the terminal session for these changes to take effect. \n

echo The GitHub environment can be activated using the command...
echo -e	"\t \$ conda activate environment_name"
echo A conda virtual environment can be deactivated using the command...
echo -e	"\t \$ conda deactivate"
