#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load MultiQC; module list

for x in 'VanEspen_2021' 'Shkoporov_2019'
do

    multiqc \
        --module fastqc \
        --filename ${x}_raw_reads_MultiQC.html \
        clean_reads_${x}/*_raw_*_fastqc.zip

    multiqc \
        --module fastqc \
        --filename ${x}_clean_reads_MultiQC.html \
        clean_reads_${x}/*_kneaddata_paired_*_fastqc.zip

done

chmod 440 *_reads_MultiQC.html
chmod 440 *_reads_MultiQC_data/*
chmod 750 *_reads_MultiQC_data
