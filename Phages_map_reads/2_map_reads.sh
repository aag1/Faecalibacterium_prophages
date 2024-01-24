#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=60gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-7
#SBATCH --export=NONE
#SBATCH --get-user-env=L



if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]; then phage='Tulp'; strain='FM7_S1'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 2 ]]; then phage='Roos'; strain='FM8_S2'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 3 ]]; then phage='Pioen'; strain='HTF-238_S3'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 4 ]]; then phage='Aster'; strain='HTF-238_S3'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 5 ]]; then phage='Lelie'; strain='L2-61_S10'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 6 ]]; then phage='CP2'; strain='FM7_S1'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 7 ]]; then phage='CP11'; strain='HTF-128_S3'; fi



SAMPLES=( $(awk -F'\t' -v strain="${strain}" '{if (($2 == strain) && ($7 == "Yes")) print $1}' ../Induction_raw_reads/Novogene_summary_table.txt) )



cd MAPPED_READS

for sample_id in "${SAMPLES[@]}" "${strain}"
do

    if [[ "${sample_id}" == "${strain}" ]]
    then
        path=../../Fprau_clean_reads/
    else
        path=../../Induction_clean_reads/CLEAN_READS/
    fi



    # ------------------------------ map reads ------------------------------ #
    base=${sample_id}_mapped_to_${phage}

    module purge; module load Bowtie2 SAMtools; module list

    bowtie2 \
        -x ../BOWTIE2_INDEX/${phage}_genome \
        -1 ${path}/${sample_id}_kneaddata_paired_1.fastq.gz \
        -2 ${path}/${sample_id}_kneaddata_paired_2.fastq.gz \
        --no-unal --threads ${SLURM_CPUS_PER_TASK} |
        samtools sort -@ $((${SLURM_CPUS_PER_TASK}-1)) - > ${base}.sorted.bam

    samtools index -@ $((${SLURM_CPUS_PER_TASK}-1)) ${base}.sorted.bam



    # ------------------------------ coverage depth & breadth ------------------------------ #
    samtools flagstat ${base}.sorted.bam > ${base}.flagstat.txt

    samtools depth ${base}.sorted.bam > ${base}.depth.txt

    module purge; module load BEDTools; module list

    bedtools coverage \
        -a ../BED_FILES/${phage}_genome.bed \
        -b ${base}.sorted.bam > ${base}.coverage.txt

    chmod 440 ${base}*

done
