#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# ------------------------------ get Pfam 36.0 seed RVT_1 MSA ------------------------------ #
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.seed.gz


module purge; module load HMMER; module list

esl-afetch \
    Pfam-A.seed.gz \
    RVT_1 \
    > RVT_1_seed.stockholm

esl-reformat \
    afa \
    RVT_1_seed.stockholm \
    > RVT_1_seed.fasta


rm Pfam-A.seed.gz RVT_1_seed.stockholm




# ------------------------------ get phage RTs ------------------------------ #
if [[ -f phages_rt.fasta ]]; then rm -f phages_rt.fasta; fi


module purge; module load seqtk; module list

for p in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie' 'CP2' 'CP4' 'CP5' 'CP6' 'CP7' 'CP11'
do

    if [[ "${p}" == "CP"* ]]
    then
        dir=../CPs_proteomes
    else
        dir=$(pwd)
    fi


    awk -F'\t' '$3 == "RVT_1" {print $1}' ${dir}/${p}_annot.txt > ${p}_rt.id

    if [[ ! -s ${p}_rt.id ]]; then continue; fi


    seqtk subseq \
        -l 80 \
        ${dir}/${p}_proteome.fasta \
        ${p}_rt.id \
        >> phages_rt.fasta

done

rm *_rt.id




# ------------------------------ add phage RTs to the Pfam MSA ------------------------------ #
module purge; module load MAFFT; module list

mafft \
    --amino \
    --add phages_rt.fasta \
    --thread ${SLURM_CPUS_PER_TASK} \
    RVT_1_seed.fasta \
    > RVT_1_seed_phages.fasta


rm RVT_1_seed.fasta phages_rt.fasta
chmod 440 RVT_1_seed_phages.fasta
