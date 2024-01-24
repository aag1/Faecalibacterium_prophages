#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-254
#SBATCH --export=NONE
#SBATCH --get-user-env=L




SAMPLE_ID=$(for x in ../VLP_reads/*_VanEspen_2021.tsv; do awk -F'\t' 'NR > 1 {print $4}' ${x}; done | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo '# ------------------------------ WORKING WITH '${SAMPLE_ID}' SAMPLE ------------------------------ #'

ARR=( $(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt) )

cd mapped_VanEspen_2021_reads




for STRAIN_ID in "${ARR[@]}"
do

    base=${STRAIN_ID}_VanEspen_2021_${SAMPLE_ID}_reads



    # ------------------------------ map reads ------------------------------ #
    module purge; module load Bowtie2 SAMtools; module list

    bowtie2 \
        -x ../BOWTIE2_INDEX/${STRAIN_ID}_contigs \
        -1 ../../VLP_reads/clean_reads_VanEspen_2021/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
        -2 ../../VLP_reads/clean_reads_VanEspen_2021/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
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

done
