#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
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



len=$(cat BED_FILES/${phage}_genome.bed | cut -d$'\t' -f3)

SAMPLES=( $(awk -F'\t' -v strain="${strain}" '{if (($2 == strain) && ($7 == "Yes")) print $1}' ../Induction_raw_reads/Novogene_summary_table.txt) )



cd NT_VARIATION

if [[ "${phage}" == "CP"* ]]
then
    genome_fasta=../../CPs_refine_coo/${phage}_refined.fna
else
    genome_fasta=../../Induction_assembly/${phage}.fasta
fi

module purge; module load Anaconda3; module list
conda activate pysamstats; conda list

for sample_id in "${SAMPLES[@]}" "${strain}"
do

    pos_cov_10=$(awk -F'\t' '$3 >= 10' ../MAPPED_READS/${sample_id}_mapped_to_${phage}.depth.txt | wc -l)

    boo=$(echo "${pos_cov_10} / ${len} >= 0.95" | bc -l)

    if [[ ${boo} -eq 0 ]]; then continue; fi


    pysamstats \
        --type variation \
        --max-depth 100000 \
        --fasta ${genome_fasta} \
        ../MAPPED_READS/${sample_id}_mapped_to_${phage}.sorted.bam \
        > ${phage}_${sample_id}_nt_variation.txt


    chmod 440 ${phage}_${sample_id}_nt_variation.txt

done
