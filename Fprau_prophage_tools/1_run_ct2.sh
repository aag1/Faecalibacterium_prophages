#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=20gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L



ARR=( $(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt) )

cd ct2_results

module purge; module load Anaconda3; module list
source activate cenote-taker2_env; conda list

for STRAIN_ID in "${ARR[@]}"
do

    STRAIN_ID2=$(echo ${STRAIN_ID} | sed 's/-/_/g')
    cp ../../Fprau_assemblies/${STRAIN_ID}_contigs.fasta ${STRAIN_ID2}_contigs.fasta


    python ${HOME}/SOFTWARE/Cenote-Taker2/run_cenote-taker2.py \
        --contigs ${STRAIN_ID2}_contigs.fasta \
        --run_title ${STRAIN_ID2} \
        --prune_prophage True \
        --virus_domain_db virion \
        --lin_minimum_hallmark_genes 1 \
        --mem 20 \
        --cpu ${SLURM_CPUS_PER_TASK}


    for f in 'final_combined_virus_sequences_'${STRAIN_ID2}'.fna' ${STRAIN_ID2}'_CONTIG_SUMMARY.tsv' ${STRAIN_ID2}'_PRUNING_INFO_TABLE.tsv'
    do
        mv ${STRAIN_ID2}/${f} .
        chmod 440 ${f}
    done

    rm -rf dummy_template.sbt ${STRAIN_ID2}_contigs* ${STRAIN_ID2}/

done

cd ../



module purge; module load R; module list

Rscript summarize_ct2.R

chmod 440 Fprau_prophages_ct2.txt
