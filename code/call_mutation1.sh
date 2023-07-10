#!/bin/bash

#SBATCH --job-name r     # Name for your job
#SBATCH --output r.out    # Standard out goes to this file
#SBATCH --error r.err     # Standard err goes to this file
#SBATCH -p short
#SBATCH --cpus-per-task=8

module load R/4.0.3-foss-2020b
Rscript code/calling_mRNA1.r






















