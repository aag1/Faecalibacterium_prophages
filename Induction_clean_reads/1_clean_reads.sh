#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-23
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# ------------------------------ sample ID ------------------------------ #
SAMPLE_ID=$(awk -F'\t' '{if ($7 == "Yes") print $1}' ../Induction_raw_reads/Novogene_summary_table.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo '# ------------------------------ WORKING WITH '${SAMPLE_ID}' SAMPLE ------------------------------ #'

cd CLEAN_READS

DD='../../Induction_raw_reads/X204SC23012117-Z01-F004/01.RawData/'${SAMPLE_ID}/
F1=$(ls ${DD}/${SAMPLE_ID}_*_1.fq.gz)
F2=$(ls ${DD}/${SAMPLE_ID}_*_2.fq.gz)
echo -e "\n\nWorking with raw read files:\n${F1}\n${F2}\n\n"




# ------------------------------ files check ------------------------------ #
WD=$(pwd)
cd ${DD}
md5sum *fq.gz > ${WD}/${SAMPLE_ID}_MD5.txt
cd ${WD}
cmp --silent ${DD}/MD5.txt ${SAMPLE_ID}_MD5.txt || echo 'Error: read files were corrupted!'




# ------------------------------ raw reads QC ------------------------------ #
module purge; module load FastQC; module list

fastqc \
    --threads 2 \
    --outdir $(pwd) \
    ${F1} \
    ${F2}

base=$(basename ${F1} _1.fq.gz)

for i in 1 2
do
    for e in 'html' 'zip'
    do
        mv ${base}_${i}_fastqc.${e} ${SAMPLE_ID}_raw_${i}_fastqc.${e}
    done
done




# ------------------------------ clean reads ------------------------------ #
module purge; module load BBMap; module list

bbduk.sh \
    in1=${F1} \
    in2=${F2} \
    out1=${SAMPLE_ID}_no_adapters_1.fq \
    out2=${SAMPLE_ID}_no_adapters_2.fq \
    ref='../../Induction_raw_reads/NEBNext_adapters.fa' \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    threads=${SLURM_CPUS_PER_TASK} \
    -Xmx20g



module purge; module load Anaconda3; module list
source activate KneadData; conda list

trimmomatic_path=${HOME}/SOFTWARE/Trimmomatic-0.33

kneaddata \
    --input ${SAMPLE_ID}_no_adapters_1.fq \
    --input ${SAMPLE_ID}_no_adapters_2.fq \
    --reference-db ../../human_genome_ref/ \
    --output-prefix ${SAMPLE_ID}_kneaddata \
    --output ${SAMPLE_ID}_kneaddata_out \
    --trimmomatic ${trimmomatic_path} \
    --trimmomatic-options "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
    --sequencer-source none \
    --bypass-trf \
    --reorder \
    --threads ${SLURM_CPUS_PER_TASK}

conda deactivate



rm ${SAMPLE_ID}_no_adapters_*.fq
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
