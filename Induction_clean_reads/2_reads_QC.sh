#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load MultiQC; module list

multiqc \
    --module fastqc \
    --filename Raw_reads_MultiQC.html \
    CLEAN_READS/*_raw_*_fastqc.zip

multiqc \
    --module fastqc \
    --filename Clean_reads_MultiQC.html \
    CLEAN_READS/*_kneaddata_paired_*_fastqc.zip



chmod 440 *_reads_MultiQC.html
chmod 440 *_reads_MultiQC_data/*
chmod 750 *_reads_MultiQC_data
