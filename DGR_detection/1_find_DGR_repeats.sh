#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




ARR=('Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie' 'CP2' 'CP4' 'CP5' 'CP6' 'CP7' 'CP11')




module purge; module load BLAST+; module list

if [[ -d BLASTN_out ]]; then rm -rf BLASTN_out; fi
mkdir BLASTN_out; chmod 750 BLASTN_out; cd BLASTN_out

for p in "${ARR[@]}"
do

    if [[ "${p}" == "CP"* ]]
    then
        f=../../CPs_refine_coo/${p}_refined.fna
    else
        f=../../Induction_assembly/${p}.fasta
    fi


    blastn \
    -task 'blastn' \
    -evalue 0.001 \
    -query ${f} \
    -subject ${f} \
    -out ${p}_repeats.11.txt \
    -outfmt 11


    blast_formatter \
        -archive ${p}_repeats.11.txt \
        -out ${p}_repeats.0.txt \
        -outfmt 0


    blast_formatter \
        -archive ${p}_repeats.11.txt \
        -out ${p}_repeats.6.txt \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'


    blast_formatter \
        -archive ${p}_repeats.11.txt \
        -out ${p}_repeats.6a.txt \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident qseq sseq'

done

cd ../




module purge; module load R/4.2.2-foss-2022b; module list

for p in "${ARR[@]}"
do

    Rscript find_DGR_repeats.R ${p}

done




chmod 440 BLASTN_out/* *_DGR_repeats.txt
