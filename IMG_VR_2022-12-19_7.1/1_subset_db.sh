#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



base='IMGVR_above10kb_DTR_noRefSeq'



awk \
    -F'\t' \
    '(NR > 1) && ($7 > 10000) && ($8 == "Direct terminal repeat") && ($19 !~ /^RefSeq/) {if ($4 == "whole") {print $1"|"$2"|"$3} else {print $1"|"$2"|"$3"|"$4}}' \
    IMGVR_all_Sequence_information-high_confidence.tsv > ${base}.ids



module purge; module load seqtk; module list

seqtk subseq \
    -l 80 \
    IMGVR_all_nucleotides-high_confidence.fna.gz \
    ${base}.ids > ${base}.fasta



chmod 440 ${base}*
