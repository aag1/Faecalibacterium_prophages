#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=60gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load BWA SAMtools; module list

bwa index \
    ../Fprau_assemblies/FM7_S1_contigs.fasta \
    -p $(pwd)/FM7_S1_contigs



for SAMPLE_ID in 'PhIFM7GN02A9' 'PhIFM7_102A3' 'PhIFM70202A4'
do

    base=${SAMPLE_ID}_mapped_to_FM7_S1_bwa

    bwa mem \
        FM7_S1_contigs \
        ../Induction_clean_reads/CLEAN_READS/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
        ../Induction_clean_reads/CLEAN_READS/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
        -t ${SLURM_CPUS_PER_TASK} |
    samtools sort -@ $((${SLURM_CPUS_PER_TASK}-1)) - > ${base}.sorted.bam

    samtools index -@ $((${SLURM_CPUS_PER_TASK}-1)) ${base}.sorted.bam

    chmod 440 ${base}*

done
