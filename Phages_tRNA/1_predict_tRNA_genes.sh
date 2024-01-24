#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load Anaconda3; module list
conda activate tRNAscan; conda list

for p in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie' 'CP2' 'CP4' 'CP5' 'CP6' 'CP7' 'CP11'
do

    if [[ "${p}" == "CP"* ]]
    then
        f=../CPs_refine_coo/${p}_refined.fna
    else
        f=../Induction_assembly/${p}.fasta
    fi


    o=${p}_tRNA_genes.txt


    tRNAscan-SE \
        -B \
        --output ${o} \
        --thread ${SLURM_CPUS_PER_TASK} \
        --forceow \
        ${f}


    if [[ ! -s ${o} ]]; then rm ${o}; fi

done

chmod 440 *_tRNA_genes.txt
