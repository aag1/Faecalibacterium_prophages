#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-455%100
#SBATCH --export=NONE
#SBATCH --get-user-env=L




SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../../IBD_reads/selected_ibd_samples.txt)

echo '# ------------------------------ WORKING WITH '${SAMPLE_ID}' SAMPLE ------------------------------ #'

cd IBD_mapped_reads




module purge; module load Bowtie2 SAMtools; module list

bowtie2 \
    -x ../five_phages \
    -1 ../../../IBD_reads/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -2 ../../../IBD_reads/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
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
