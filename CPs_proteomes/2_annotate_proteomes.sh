#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-15
#SBATCH --export=NONE
#SBATCH --get-user-env=L



ID='CP'${SLURM_ARRAY_TASK_ID}



module purge; module load HMMER; module list

hmmsearch \
    --max \
    -E 0.001 \
    --cpu ${SLURM_CPUS_PER_TASK} \
    --domtblout ${ID}_vs_Pfam.txt \
    -o ${ID}_vs_Pfam.out \
    ../Pfam_36.0/Pfam-A.hmm \
    ${ID}_proteome.fasta



module purge; module load R/4.2.2-foss-2022b; module list

Rscript parse_hmmer.R ${ID} ../Pfam_36.0/PfamA_36.0_domains.txt



chmod 440 ${ID}_vs_Pfam.* ${ID}_annot.txt
