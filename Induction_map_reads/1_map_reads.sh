#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A_%a.out
#SBATCH --mem=60gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-23
#SBATCH --export=NONE
#SBATCH --get-user-env=L




SAMPLE_ID=$(awk -F'\t' '{if ($7 == "Yes") print $1}' ../Induction_raw_reads/Novogene_summary_table.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")

STRAIN_ID=$(awk -F'\t' '{if ($7 == "Yes") print $2}' ../Induction_raw_reads/Novogene_summary_table.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")

base=${SAMPLE_ID}_mapped_to_${STRAIN_ID}

cd MAPPED_READS




# ------------------------------ map reads ------------------------------ #
module purge; module load Bowtie2 SAMtools; module list

bowtie2 \
    -x ../../Fprau_map_reads/BOWTIE2_INDEX/${STRAIN_ID}_contigs \
    -1 ../../Induction_clean_reads/CLEAN_READS/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -2 ../../Induction_clean_reads/CLEAN_READS/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    --no-unal --threads ${SLURM_CPUS_PER_TASK} |
samtools sort -@ $((${SLURM_CPUS_PER_TASK}-1)) - > ${base}.sorted.bam

samtools index -@ $((${SLURM_CPUS_PER_TASK}-1)) ${base}.sorted.bam




# ------------------------------ coverage depth & breadth ------------------------------ #
samtools flagstat ${base}.sorted.bam > ${base}.flagstat.txt

samtools depth ${base}.sorted.bam > ${base}.depth.txt

module purge; module load BEDTools; module list

bedtools coverage \
    -a ../../Fprau_map_reads/BED_FILES/${STRAIN_ID}_contigs.bed \
    -b ${base}.sorted.bam > ${base}.coverage.txt

chmod 440 ${base}*
