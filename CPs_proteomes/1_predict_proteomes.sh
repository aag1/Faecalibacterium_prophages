#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load Anaconda3; module list
conda activate mgcod; conda list

for i in {1..15}
do

    mgcod.py \
        --path_to_genome ../CPs_refine_coo/CP${i}_refined.fna \
        --isoforms \
        --path_to_mgm_predictions CP${i}_mgcod \
        --path_to_output CP${i}_mgcod \
        --logfile CP${i}_mgcod/CP${i}_mgcod.log \
        --amino_acids \
        --delete \
        --verbose


    awk \
        -F'\t' \
        -v OFS='\t' \
        '($1 !~ /^#/) && ($1 !~ /^$/) && ($6 != ".") {print $9,$4,$5,$7}' \
        CP${i}_mgcod/CP${i}_refined.gff |
        sed -E "s/^gene_id ([0-9\.]+);[^\t]+(\t.+)$/CP${i}_\1\2/" \
        > CP${i}_proteome.txt


    sed \
        -E "s/^>.+ gene_id ([0-9\.]+); .+$/>CP${i}_\1/" \
        CP${i}_mgcod/proteins_aa_CP${i}_refined.fasta \
        > CP${i}_proteome.fasta

done



grep 'genetic code ' *_mgcod/*log



mkdir Mgcod_OUT; chmod 750 Mgcod_OUT; mv CP*_mgcod Mgcod_OUT
chmod 440 *_proteome.txt *_proteome.fasta
