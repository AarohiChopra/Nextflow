#!/bin/bash
#SBATCH --job-name=Nextflow         	# Job name
#SBATCH --nodes=1                   		# Number of nodes
#SBATCH --ntasks=1                  		# Number of tasks (processes)
#SBATCH --cpus-per-task=5           		# Number of CPU cores per task
#SBATCH --mem=10G                   		# Total memory for the job
#SBATCH --time=01:00:00             		# Time limit hrs:min:sec
#SBATCH --output=/home/achopra/BPA_Alt_Human/BPA/Nextflow/scratch/Nextflow_log%j.out  		# Standard output and error log

# README sbatch runNextflow.sh

source activate /home/achopra/miniconda3/envs/nextFlow  # Load Nextflow module, if available on your cluster

nextflow run pipeline2.nf -resume
