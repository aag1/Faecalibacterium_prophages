#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



while read -r seq1 seq2 seq2_strand seq1_len seq2_len seq2_db
do

    if [[ "${seq1}" == "seq1" ]]; then continue; fi

    SEQ2=$(echo ${seq2} | cut -d'|' -f1)

    pair_id=${seq1}'___'${SEQ2}


    ### first sequence
    cat ${seq1}_raw.fna > ${pair_id}.fasta


    ### second sequence
    if [[ "${seq2_db}" == "viral_refseq" ]]
    then
        db_fasta=../Viral_RefSeq_220/viral_refseq_complete_above10kb_noRiboviria.fasta
    else
        db_fasta=../IMG_VR_2022-12-19_7.1/IMGVR_above10kb_DTR_noRefSeq.fasta
    fi


    module purge; module load seqtk; module list
    echo ${seq2} > seq2.id
    seqtk subseq -l 80 ${db_fasta} seq2.id > seq2.fasta
    sed -i -E 's/^(>[^\|]+)\|.+$/\1/' seq2.fasta


    module purge; module load EMBOSS; module list
    if [[ "${seq2_strand}" == "r" ]]; then revseq seq2.fasta seq2.fasta; fi


    cat seq2.fasta >> ${pair_id}.fasta


    ### dotplot
    polydot ${pair_id}.fasta -wordsize 12 -graph pdf -dumpfeat -outfeat ${pair_id}.gff
    mv polydot.pdf ${pair_id}.pdf

done < pairs_for_dotplots.txt



rm seq2.id seq2.fasta *___*fasta
chmod 440 *___*
mkdir dotplots_data; chmod 750 dotplots_data
mv *___* dotplots_data/



module purge; module load R; module list
Rscript dotplots.R
chmod 440 dotplots.pdf
