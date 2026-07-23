#!/bin/bash

#Request cores. In this case - 2 cores will be requested
#SBATCH --cpus-per-task=1

#Request memory per thread. Here it's 1 GB per requested core
#SBATCH --mem-per-cpu=30G

# Set a job name
#SBATCH --job-name=CODEX_coloc

#Specify time before the job is killed by scheduler (in case is hangs). In this case - 24 hours
#SBATCH --time=24:00:00 

#Declare the merged STDOUT/STDERR file
#SBATCH --output=output.%J.out

#Run something
#conda activate new_NicheSphere_env
python ./019_CODEX_coloc_radius.py
