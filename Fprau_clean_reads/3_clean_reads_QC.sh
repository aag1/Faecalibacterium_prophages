#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load FastQC; module list

ARR=( $(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt) )

for x in "${ARR[@]}"
do
    fastqc \
        --threads 2 \
        --outdir $(pwd) \
        ${x}_kneaddata_paired_1.fastq.gz \
        ${x}_kneaddata_paired_2.fastq.gz
done



module purge; module load MultiQC; module list

multiqc \
    --module fastqc \
    --filename Fprau_clean_reads_MultiQC.html \
    *_fastqc.zip



mkdir Fprau_clean_reads_FastQC
chmod 750 Fprau_clean_reads_FastQC
mv *_fastqc* Fprau_clean_reads_FastQC/
