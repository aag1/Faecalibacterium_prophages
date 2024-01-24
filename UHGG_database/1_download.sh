#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=0-14
#SBATCH --export=NONE
#SBATCH --get-user-env=L



# wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/README_v2.0.2.txt
# wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/genomes-all_metadata.tsv



ARR=( $(awk -F'\t' '((NR > 1) && (($1 == $14) || ($2 == "Isolate"))) {print $20}' genomes-all_metadata.tsv) )



from=$((${SLURM_ARRAY_TASK_ID} * 1000))

arr=("${ARR[@]:${from}:1000}")



cd fna_files

for f in "${arr[@]}"
do

    g=$(basename ${f} .gff.gz)

    echo 'Working with '${g}' genome...'

    wget --no-verbose ${f}

    gunzip -c ${g}.gff.gz | sed -n '/^##FASTA$/,$p' | sed '1d' > ${g}.fna

    rm ${g}.gff.gz

    chmod 440 ${g}.fna

done
