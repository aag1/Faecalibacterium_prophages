#!/bin/bash
#SBATCH --job-name=job5
#SBATCH --output=job5_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



mkdir Blast_OUT; chmod 750 Blast_OUT; cd Blast_OUT

for i in 2 3 4 10 11 12
do

    phage=$(sed -n "${i}p" ../induced_phages_2.txt | cut -d$'\t' -f1)

    contig=$(sed -n "${i}p" ../induced_phages_2.txt | cut -d$'\t' -f4)

    sample=$(sed -n "${i}p" ../induced_phages_2.txt | cut -d$'\t' -f7)

    strain=$(sed -n "${i}p" ../induced_phages_2.txt | cut -d$'\t' -f8)



    module purge; module load seqtk; module list

    echo ${contig} > contig.id

    seqtk subseq \
        -l 80 \
        ../induced_phages.fasta \
        contig.id > contig.fasta



    module purge; module load BLAST+; module list

    makeblastdb \
        -in ../../Fprau_assemblies/${strain}_contigs.fasta \
        -dbtype nucl \
        -out ${strain}

    blastn \
        -task 'blastn' \
        -db ${strain} \
        -query contig.fasta \
        -evalue 0.001 \
        -num_threads ${SLURM_CPUS_PER_TASK} \
        -out ${phage}'_'${sample}'_vs_'${strain}'.txt' \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'



    chmod 440 ${phage}'_'${sample}'_vs_'${strain}'.txt'
    rm contig.id contig.fasta ${strain}.n*

done
