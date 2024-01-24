#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L




if [[ -f five_phages.fasta ]]; then rm -f five_phages.fasta; fi

if [[ -f five_phages.bed ]]; then rm -f five_phages.bed; fi

for p in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie'
do

    cat ../Induction_assembly/${p}.fasta >> five_phages.fasta

    cat ../Phages_map_reads/BED_FILES/${p}_genome.bed >> five_phages.bed

done




module purge; module load Bowtie2; module list

bowtie2-build \
    five_phages.fasta \
    five_phages \
    --threads ${SLURM_CPUS_PER_TASK}

chmod 440 five_phages*
