#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load prodigal; module list

ARR=( $(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt) )

for STRAIN_ID in "${ARR[@]}"
do
    prodigal \
        -p single \
        -i ../Fprau_assemblies/${STRAIN_ID}_contigs.fasta \
        -f gff -o ${STRAIN_ID}_proteome.gff
done

chmod 440 *_proteome.gff
