#!/bin/bash
#SBATCH --job-name=job3.1
#SBATCH --output=job3.1_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-23
#SBATCH --export=NONE
#SBATCH --get-user-env=L



STRAIN_ID=$(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo '# ------------------------------ WORKING WITH '${STRAIN_ID}' STRAIN ------------------------------ #'

cd geNomad_results



module purge; module load Anaconda3; module list
source activate geNomad; conda list

genomad end-to-end \
    --cleanup \
    --threads ${SLURM_CPUS_PER_TASK} \
    ../../Fprau_assemblies/${STRAIN_ID}_contigs.fasta \
    ${STRAIN_ID} \
    ${HOME}/SOFTWARE/genomad_db
