#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load SAMtools; module list

samtools view PhIFM7_102A3_mapped_to_FM7_S1_bwa.sorted.bam |
awk -F'\t' '$17 ~ /^SA:Z:NODE_11_/' > PhIFM7_102A3_split_info.txt



module purge; module load R; module list

Rscript get_split_coo.R

Rscript plot_transpos_sites.R



chmod 440 PhIFM7_102A3_split_* Tulp_transpos_sites.pdf
