#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-5
#SBATCH --export=NONE
#SBATCH --get-user-env=L



if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]; then p='Tulp'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 2 ]]; then p='Roos'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 3 ]]; then p='Pioen'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 4 ]]; then p='Aster'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 5 ]]; then p='Lelie'; fi



module purge; module load HMMER; module list

hmmsearch \
    --max \
    -E 0.001 \
    --cpu ${SLURM_CPUS_PER_TASK} \
    --domtblout ${p}_vs_Pfam.txt \
    -o ${p}_vs_Pfam.out \
    ../Pfam_36.0/Pfam-A.hmm \
    ${p}_proteome.fasta



module purge; module load R/4.2.2-foss-2022b; module list

Rscript ../CPs_proteomes/parse_hmmer.R ${p} ../Pfam_36.0/PfamA_36.0_domains.txt



chmod 440 ${p}_vs_Pfam* ${p}_annot.txt
