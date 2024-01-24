#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-6
#SBATCH --export=NONE
#SBATCH --get-user-env=L



if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]; then p='Tulp'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 2 ]]; then p='Roos'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 3 ]]; then p='Pioen'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 4 ]]; then p='Lelie'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 5 ]]; then p='CP2'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 6 ]]; then p='CP11'; fi



if [[ "${p}" == "CP"* ]]
then
    f=../CPs_refine_coo/${p}_refined.fna
else
    f=../Induction_assembly/${p}.fasta
fi



module purge; module load BLAST+; module list

blastn \
    -task 'blastn' \
    -db ../Viral_RefSeq_220/viral_refseq_complete_above10kb_noRiboviria \
    -query ${f} \
    -evalue 0.001 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -out ${p}_vs_viral_refseq.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'



blastn \
    -task 'blastn' \
    -db ../IMG_VR_2022-12-19_7.1/IMGVR_above10kb_DTR_noRefSeq \
    -query ${f} \
    -evalue 0.001 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -out ${p}_vs_IMGVR.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'



N=$(awk -F'\t' '$10 >= 95 {print $2}' ${p}_vs_*.txt | sort | uniq | wc -l)
echo ${N}' cognate genomes found for '${p}'!'



chmod 440 ${p}_vs_*.txt
