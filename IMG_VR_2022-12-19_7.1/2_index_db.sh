#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



base='IMGVR_above10kb_DTR_noRefSeq'



module purge; module load BLAST+; module list

makeblastdb \
    -in ${base}.fasta \
    -dbtype nucl \
    -out ${base}
