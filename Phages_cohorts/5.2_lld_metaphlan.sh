#!/bin/bash
#SBATCH --job-name=job5.2
#SBATCH --output=job5.2_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-135
#SBATCH --export=NONE
#SBATCH --get-user-env=L




N=$((${SLURM_ARRAY_TASK_ID} + 1000))

SAMPLE_ID=$(cut -d$'\t' -f1 ../../LLD_reads/LLD_raw_reads_number_sele.txt | sed "${N}q;d")

echo '# ------------------------------ WORKING WITH '${SAMPLE_ID}' SAMPLE ------------------------------ #'

cd LLD_metaphlan




module purge; module load Anaconda3; module list
source activate MetaPhlAn_4.1.1; conda list

metaphlan \
    --input_type fastq \
    ../../../LLD_reads/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz,../../../LLD_reads/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    --bowtie2db ../../../MetaPhlAn_4.1.1_DB \
    --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
    --bowtie2out ${SAMPLE_ID}_metaphlan.bz2 \
    --output_file ${SAMPLE_ID}_metaphlan.txt \
    --nproc ${SLURM_CPUS_PER_TASK}

rm ${SAMPLE_ID}_metaphlan.bz2
chmod 440 ${SAMPLE_ID}_metaphlan.txt




echo '# ------------------------------ FINISHED WORKING! ------------------------------ #'
