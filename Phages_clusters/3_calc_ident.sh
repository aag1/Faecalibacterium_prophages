#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



if [[ -f pairwise_identity.txt ]]; then rm -f pairwise_identity.txt; fi

while read -r seq1 seq2
do

    ### first sequence
    if [[ -f ../Induction_assembly/${seq1}.fasta ]];
    then
        f1=../Induction_assembly/${seq1}.fasta
    else
        f1=../CPs_refine_coo/${seq1}_refined.fna
    fi


    ### second sequence
    module purge; module load seqtk; module list

    echo ${seq2} > seq2.id

    seqtk subseq -l 80 ../Viral_RefSeq_220/viral.1.1.genomic.fna seq2.id > seq2.fasta



    ### prepare for alignment
    if [[ "${seq1}" == "Tulp" ]]
    then

        cat ${f1} > ${seq1}_${seq2}_for_ali.fasta

        cat seq2.fasta >> ${seq1}_${seq2}_for_ali.fasta

    else

        ### first sequence 5'-end
        echo -e ${seq1}'\t0\t500' > seq1_5p.bed

        seqtk subseq -l 80 ${f1} seq1_5p.bed > seq1_5p.fasta



        ### find new 5'-end
        module purge; module load BLAST+; module list

        blastn \
            -task 'blastn' \
            -query seq1_5p.fasta \
            -subject seq2.fasta \
            -out ${seq1}_5p_vs_${seq2}.txt \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'



        ### apply new 5'-end
        module purge; module load R/4.2.2-foss-2022b; module list

        Rscript new_5prime.R \
            ${f1} \
            seq2.fasta \
            ${seq1}_5p_vs_${seq2}.txt \
            ${seq1}_${seq2}_for_ali.fasta

    fi



    ### align genomes
    module purge; module load MAFFT; module list

    mafft \
        --nuc \
        --maxiterate 1000 \
        --thread ${SLURM_CPUS_PER_TASK} \
        ${seq1}_${seq2}_for_ali.fasta > ${seq1}_${seq2}_ali.fasta



    ### calculate identity
    module purge; module load R/4.2.2-foss-2022b; module list

    Rscript calc_ident.R ${seq1} ${seq2} >> pairwise_identity.txt

done < pairs_for_dotplots.txt



rm seq1_5p* seq2* *_for_ali.fasta
chmod 440 *_5p_vs_*.txt *_ali.fasta pairwise_identity.txt
