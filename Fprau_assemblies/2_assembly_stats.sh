#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load BBMap; module list

statswrapper.sh *_contigs.fasta > stats_contigs.txt

statswrapper.sh *_scaffolds.fasta > stats_scaffolds.txt

statswrapper.sh ../Fprau_assemblies_Lei/*fna > stats_Lei.txt



module purge; module load R; module list

Rscript compare_stats.R



chmod 440 stats_*txt scaffolds_my_vs_Lei.pdf
