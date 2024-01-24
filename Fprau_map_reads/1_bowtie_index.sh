#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L




ARR=( $(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt) )




# ------------------------------ bowtie2 indices ------------------------------ #
module purge; module load Bowtie2; module list

mkdir BOWTIE2_INDEX; chmod 750 BOWTIE2_INDEX; cd BOWTIE2_INDEX

for STRAIN_ID in "${ARR[@]}"
do
    bowtie2-build \
        ../../Fprau_assemblies/${STRAIN_ID}_contigs.fasta \
        ${STRAIN_ID}_contigs \
        --threads ${SLURM_CPUS_PER_TASK}
done

cd ../




# ------------------------------ bed files ------------------------------ #
module purge; module load EMBOSS; module list

mkdir BED_FILES; chmod 750 BED_FILES; cd BED_FILES

for STRAIN_ID in "${ARR[@]}"
do
    infoseq \
        -auto \
        -nocolumns -noheading \
        -delimiter $'\t0\t' \
        -only -name -length \
        ../../Fprau_assemblies/${STRAIN_ID}_contigs.fasta \
        > ${STRAIN_ID}_contigs.bed
done

chmod 440 *_contigs.bed
