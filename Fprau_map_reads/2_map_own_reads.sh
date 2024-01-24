#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-23
#SBATCH --export=NONE
#SBATCH --get-user-env=L




STRAIN_ID=$(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo '# ------------------------------ WORKING WITH '${STRAIN_ID}' STRAIN ------------------------------ #'

base=${STRAIN_ID}_own_reads

cd mapped_own_reads




# ------------------------------ map reads ------------------------------ #
module purge; module load Bowtie2 SAMtools; module list

bowtie2 \
    -x ../BOWTIE2_INDEX/${STRAIN_ID}_contigs \
    -1 ../../Fprau_clean_reads/${STRAIN_ID}_kneaddata_paired_1.fastq.gz \
    -2 ../../Fprau_clean_reads/${STRAIN_ID}_kneaddata_paired_2.fastq.gz \
    --no-unal --threads ${SLURM_CPUS_PER_TASK} |
samtools sort -@ $((${SLURM_CPUS_PER_TASK}-1)) - > ${base}.sorted.bam

samtools index -@ $((${SLURM_CPUS_PER_TASK}-1)) ${base}.sorted.bam




# ------------------------------ coverage depth & breadth ------------------------------ #
samtools flagstat ${base}.sorted.bam > ${base}.flagstat.txt

samtools depth ${base}.sorted.bam > ${base}.depth.txt

module purge; module load BEDTools; module list

bedtools coverage \
    -a ../BED_FILES/${STRAIN_ID}_contigs.bed \
    -b ${base}.sorted.bam > ${base}.coverage.txt

chmod 440 ${base}*
