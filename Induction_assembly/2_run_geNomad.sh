#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=20gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L



cd geNomad_OUT

module purge; module load Anaconda3; module list
source activate geNomad; conda list

for f in ../*_contigs.fasta
do

    s=$(basename ${f} _contigs.fasta)

    genomad end-to-end \
        --cleanup \
        --threads ${SLURM_CPUS_PER_TASK} \
        ${f} \
        ${s} \
        ${HOME}/SOFTWARE/genomad_db

done
