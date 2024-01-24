#!/bin/bash
#SBATCH --job-name=job1.2
#SBATCH --output=job1.2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=6-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2
#SBATCH --export=NONE
#SBATCH --get-user-env=L


# https://www.ebi.ac.uk/ena/browser/view/PRJNA722819
# https://www.ebi.ac.uk/ena/browser/view/PRJNA723467
# 'Show Column Selection' --> add/remove columns --> Download report: TSV


if [ ${SLURM_ARRAY_TASK_ID} -eq 1 ]; then tsv_file='PRJNA722819_VanEspen_2021.tsv'; else tsv_file='PRJNA723467_VanEspen_2021.tsv'; fi
ARR=( $(awk -F'\t' 'NR > 1 {print $9}' ${tsv_file}) )

cd raw_reads_VanEspen_2021

for x in "${ARR[@]}"
do

    f1=$(echo ${x} | cut -d';' -f1)
    f2=$(echo ${x} | cut -d';' -f2)

    wget --no-verbose ${f1}
    wget --no-verbose ${f2}

done
