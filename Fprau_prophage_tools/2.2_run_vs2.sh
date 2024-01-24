#!/bin/bash
#SBATCH --job-name=job2.2
#SBATCH --output=job2.2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load R; module list

Rscript summarize_vs2.R

chmod 440 Fprau_prophages_vs2.txt
