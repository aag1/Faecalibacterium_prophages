#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



# ------------------------------ bowtie2 indices ------------------------------ #
module purge; module load Bowtie2; module list

mkdir BOWTIE2_INDEX; chmod 750 BOWTIE2_INDEX; cd BOWTIE2_INDEX

for p in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie' 'CP2' 'CP11'
do

    if [[ "${p}" == "CP"* ]]
    then
        f=../../CPs_refine_coo/${p}_refined.fna
    else
        f=../../Induction_assembly/${p}.fasta
    fi


    bowtie2-build \
        ${f} \
        ${p}_genome \
        --threads ${SLURM_CPUS_PER_TASK}

done

cd ../



# ------------------------------ bed files ------------------------------ #
module purge; module load EMBOSS; module list

mkdir BED_FILES; chmod 750 BED_FILES; cd BED_FILES

for p in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie' 'CP2' 'CP11'
do

    if [[ "${p}" == "CP"* ]]
    then
        f=../../CPs_refine_coo/${p}_refined.fna
    else
        f=../../Induction_assembly/${p}.fasta
    fi


    infoseq \
        -auto \
        -nocolumns -noheading \
        -delimiter $'\t0\t' \
        -only -name -length \
        ${f} \
        > ${p}_genome.bed

done

chmod 440 *_genome.bed
