#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=0-109
#SBATCH --export=NONE
#SBATCH --get-user-env=L



ARR=( $(awk -F'\t' '$2 == "Isolate" {print $1}' genomes-all_metadata.tsv) )



from=$((${SLURM_ARRAY_TASK_ID} * 100))

arr=("${ARR[@]:${from}:100}")



cd CRISPRidentify_out

module purge; module load Anaconda3; module list
conda activate crispr_identify_env; conda list

for g in "${arr[@]}"
do

    echo 'Working with '${g}' genome...'


    dir=${g}_CRISPRidentify


    python ${HOME}/SOFTWARE/CRISPRidentify/CRISPRidentify.py \
        --file ../fna_files/${g}.fna \
        --result_folder ${dir} \
        --fasta_report True \
        --cpu ${SLURM_CPUS_PER_TASK}

    
    if [[ -s ${dir}/Complete_spacer_dataset.fasta ]]
    then
        mv ${dir}/Complete_spacer_dataset.fasta ${g}_Complete_spacer_dataset.fasta
        mv ${dir}/Complete_summary.csv ${g}_Complete_summary.csv
        chmod 440 ${g}_Complete_spacer_dataset.fasta ${g}_Complete_summary.csv
    fi


    rm -rf ${dir}

done
