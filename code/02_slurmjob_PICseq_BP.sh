#!/bin/bash

#Request cores. In this case - 25 cores will be requested
#SBATCH --cpus-per-task=25

#Request memory per thread. Here it's 30 GB per requested core
#SBATCH --mem-per-cpu=30G

# Set a job name
#SBATCH --job-name=BP4PICseq

#Specify time before the job is killed by scheduler (in case is hangs). In this case - 60 hours
#SBATCH --time=60:00:00 

#Declare the merged STDOUT/STDERR file
#SBATCH --output=output.%J.out

#Run something
#export PATH="/home/hu367653/miniconda3/bin:$PATH"

#conda activate sceasy_env
Rscript 06_PICseq_BP.R
