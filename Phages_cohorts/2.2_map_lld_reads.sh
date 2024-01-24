#!/bin/bash
#SBATCH --job-name=job2.2
#SBATCH --output=job2.2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-135
#SBATCH --export=NONE
#SBATCH --get-user-env=L




N=$((${SLURM_ARRAY_TASK_ID} + 1000))

SAMPLE_ID=$(cut -d$'\t' -f1 ../../LLD_reads/LLD_raw_reads_number_sele.txt | sed "${N}q;d")

echo '# ------------------------------ WORKING WITH '${SAMPLE_ID}' SAMPLE ------------------------------ #'

cd LLD_mapped_reads




module purge; module load Bowtie2 SAMtools; module list

bowtie2 \
    -x ../five_phages \
    -1 ../../../LLD_reads/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -2 ../../../LLD_reads/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    --no-unal --threads ${SLURM_CPUS_PER_TASK} |
    samtools sort -@ $((${SLURM_CPUS_PER_TASK}-1)) - > ${SAMPLE_ID}_sorted.bam

samtools index -@ $((${SLURM_CPUS_PER_TASK}-1)) ${SAMPLE_ID}_sorted.bam

samtools flagstat ${SAMPLE_ID}_sorted.bam > ${SAMPLE_ID}_flagstat.txt




module purge; module load BEDTools; module list

bedtools coverage \
    -a ../five_phages.bed \
    -b ${SAMPLE_ID}_sorted.bam > ${SAMPLE_ID}_coverage.txt

chmod 440 ${SAMPLE_ID}_*




echo '# ------------------------------ FINISHED WORKING! ------------------------------ #'
