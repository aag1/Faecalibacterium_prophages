#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A_%a.out
#SBATCH --mem=60gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-17
#SBATCH --export=NONE
#SBATCH --get-user-env=L



SAMPLE_ID=$(awk -F'\t' '($3 != "bacterial DNA") && ($7 == "Yes") {print $1}' ../Induction_raw_reads/Novogene_summary_table.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")



module purge; module load SPAdes; module list

spades.py \
    --metaviral \
    -1 ../Induction_clean_reads/CLEAN_READS/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz \
    -2 ../Induction_clean_reads/CLEAN_READS/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    -o ${SAMPLE_ID}_spades_out \
    --threads ${SLURM_CPUS_PER_TASK}



for x in 'contigs.fasta' 'spades.log'
do
    file=${SAMPLE_ID}_spades_out/${x}
    if [[ -f ${file} ]]; then mv ${file} ${SAMPLE_ID}_${x}; fi
done

rm -rf ${SAMPLE_ID}_spades_out
chmod 440 ${SAMPLE_ID}_*
