#!/bin/bash
#SBATCH --job-name=job3.2
#SBATCH --output=job3.2_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-254%50
#SBATCH --export=NONE
#SBATCH --get-user-env=L



SAMPLE_ID=$(for x in *_VanEspen_2021.tsv; do awk -F'\t' 'NR > 1 {print $4}' ${x}; done | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo '# ------------------------------ WORKING WITH '${SAMPLE_ID}' SAMPLE ------------------------------ #'

cd clean_reads_VanEspen_2021



# ------------------------------ raw reads QC ------------------------------ #
module purge; module load FastQC; module list

fastqc \
    --threads 2 \
    --outdir $(pwd) \
    ../raw_reads_VanEspen_2021/${SAMPLE_ID}_1.fastq.gz \
    ../raw_reads_VanEspen_2021/${SAMPLE_ID}_2.fastq.gz

for i in 1 2
do
    for e in 'html' 'zip'
    do
        mv ${SAMPLE_ID}_${i}_fastqc.${e} ${SAMPLE_ID}_raw_${i}_fastqc.${e}
    done
done



# ------------------------------ clean reads ------------------------------ #
module purge; module load Anaconda3; module list
source activate KneadData; conda list

trimmomatic_path=${HOME}/SOFTWARE/Trimmomatic-0.33

kneaddata \
    --input ../raw_reads_VanEspen_2021/${SAMPLE_ID}_1.fastq.gz \
    --input ../raw_reads_VanEspen_2021/${SAMPLE_ID}_2.fastq.gz \
    --reference-db ../../human_genome_ref/ \
    --output-prefix ${SAMPLE_ID}_kneaddata \
    --output ${SAMPLE_ID}_kneaddata_out \
    --trimmomatic ${trimmomatic_path} \
    --trimmomatic-options "ILLUMINACLIP:${trimmomatic_path}/adapters/NexteraPE-PE.fa:2:30:10:1:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
    --sequencer-source NexteraPE \
    --bypass-trf \
    --reorder \
    --threads ${SLURM_CPUS_PER_TASK}

conda deactivate

mv ${SAMPLE_ID}_kneaddata_out/${SAMPLE_ID}_kneaddata_paired_1.fastq .
mv ${SAMPLE_ID}_kneaddata_out/${SAMPLE_ID}_kneaddata_paired_2.fastq .
mv ${SAMPLE_ID}_kneaddata_out/${SAMPLE_ID}_kneaddata.log .
rm -rf ${SAMPLE_ID}_kneaddata_out



# ------------------------------ clean reads QC ------------------------------ #
module purge; module load FastQC; module list

fastqc \
    --threads 2 \
    --outdir $(pwd) \
    ${SAMPLE_ID}_kneaddata_paired_1.fastq \
    ${SAMPLE_ID}_kneaddata_paired_2.fastq

gzip ${SAMPLE_ID}_kneaddata_paired_1.fastq
gzip ${SAMPLE_ID}_kneaddata_paired_2.fastq
