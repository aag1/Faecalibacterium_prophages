#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



while read -r seq1 seq2
do

    pair_id=${seq1}'___'${seq2}


    ### first sequence
    if [[ -f ../Induction_assembly/${seq1}.fasta ]];
    then
        cat ../Induction_assembly/${seq1}.fasta > ${pair_id}.fasta
    else
        cat ../CPs_refine_coo/${seq1}_refined.fna > ${pair_id}.fasta
    fi


    ### second sequence
    module purge; module load seqtk; module list

    db_fasta=../Viral_RefSeq_220/viral.1.1.genomic.fna

    echo ${seq2} > seq2.id

    seqtk subseq -l 80 ${db_fasta} seq2.id >> ${pair_id}.fasta


    ### dotplot
    module purge; module load EMBOSS; module list

    polydot ${pair_id}.fasta -wordsize 12 -graph pdf -dumpfeat -outfeat ${pair_id}.gff

    mv polydot.pdf ${pair_id}.pdf

done < pairs_for_dotplots.txt



rm seq2.id *___*fasta
chmod 440 *___*
mkdir dotplots_data; chmod 750 dotplots_data
mv *___* dotplots_data/



module purge; module load R/4.2.2-foss-2022b; module list
Rscript dotplots.R
chmod 440 dotplots_VCs.pdf
