#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



base='CPs_induced_phages'

cat ../CPs_clusters/CPs15.fna > ${base}.fna

cat induced_phages.fasta >> ${base}.fna



module purge; module load BLAST+; module list

makeblastdb \
    -in ${base}.fna \
    -dbtype nucl \
    -out ${base}

blastn \
    -query ${base}.fna \
    -db ${base} \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 10000 \
    -out ${base}_blast.tsv \
    -num_threads ${SLURM_CPUS_PER_TASK}



module purge; module load CheckV; module list

python ../CPs_clusters/anicalc.py \
    -i ${base}_blast.tsv \
    -o ${base}_ani.tsv

python ../CPs_clusters/aniclust.py \
    --fna ${base}.fna \
    --ani ${base}_ani.tsv \
    --out ${base}_clusters.tsv \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0



rm ${base}.fna ${base}.n*
chmod 440 ${base}*
