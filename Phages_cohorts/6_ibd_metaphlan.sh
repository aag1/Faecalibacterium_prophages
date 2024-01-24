#!/bin/bash
#SBATCH --job-name=job6
#SBATCH --output=job6_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-455
#SBATCH --export=NONE
#SBATCH --get-user-env=L




SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../../IBD_reads/selected_ibd_samples.txt)

echo '# ------------------------------ WORKING WITH '${SAMPLE_ID}' SAMPLE ------------------------------ #'

cd IBD_metaphlan




module purge; module load Anaconda3; module list
source activate MetaPhlAn4; conda list

metaphlan \
    --input_type fastq \
    ../../../IBD_reads/${SAMPLE_ID}_kneaddata_paired_1.fastq.gz,../../../IBD_reads/${SAMPLE_ID}_kneaddata_paired_2.fastq.gz \
    --bowtie2db ../../../MetaPhlAn4_DB \
    --index mpa_vOct22_CHOCOPhlAnSGB_202212 \
    --bowtie2out ${SAMPLE_ID}_metaphlan.bz2 \
    --output_file ${SAMPLE_ID}_metaphlan.txt \
    --nproc ${SLURM_CPUS_PER_TASK}

rm ${SAMPLE_ID}_metaphlan.bz2
chmod 440 ${SAMPLE_ID}_metaphlan.txt




echo '# ------------------------------ FINISHED WORKING! ------------------------------ #'
