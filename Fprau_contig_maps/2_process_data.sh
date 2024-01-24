#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-23
#SBATCH --export=NONE
#SBATCH --get-user-env=L



STRAIN_ID=$(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo '# ------------------------------ WORKING WITH '${STRAIN_ID}' STRAIN ------------------------------ #'



module purge; module load R; module list

Rscript summarize_own_reads_coverage.R ${STRAIN_ID}

Rscript summarize_vlp_coverage.R ${STRAIN_ID}

Rscript calculate_nt_content.R ${STRAIN_ID}



chmod 440 ${STRAIN_ID}*RData
