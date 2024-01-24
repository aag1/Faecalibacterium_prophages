#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-15
#SBATCH --export=NONE
#SBATCH --get-user-env=L



Q='CP'${SLURM_ARRAY_TASK_ID}'_raw'



module purge; module load BLAST+; module list

blastn \
    -task 'blastn' \
    -db ../Viral_RefSeq_220/viral_refseq_complete_above10kb_noRiboviria \
    -query ${Q}.fna \
    -evalue 0.001 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -out ${Q}_vs_viral_refseq.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'



blastn \
    -task 'blastn' \
    -db ../IMG_VR_2022-12-19_7.1/IMGVR_above10kb_DTR_noRefSeq \
    -query ${Q}.fna \
    -evalue 0.001 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -out ${Q}_vs_IMGVR.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'



chmod 440 ${Q}_vs_*.txt
