#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



base='CPs15'

cat ../CPs_refine_coo/*_refined.fna > ${base}.fna



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

python anicalc.py \
    -i ${base}_blast.tsv \
    -o ${base}_ani.tsv

python aniclust.py \
    --fna ${base}.fna \
    --ani ${base}_ani.tsv \
    --out ${base}_clusters.tsv \
    --min_ani 95 \
    --min_tcov 85 \
    --min_qcov 0



rm CPs15.n*
chmod 440 CPs15*
