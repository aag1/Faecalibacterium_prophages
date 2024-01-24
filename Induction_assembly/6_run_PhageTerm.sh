#!/bin/bash
#SBATCH --job-name=job6
#SBATCH --output=job6_%A.out
#SBATCH --mem=20gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L



mkdir PhageTerm_OUT; chmod 750 PhageTerm_OUT; cd PhageTerm_OUT

module purge; module load seqtk; module list

for phage in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie'
do

    ### get phage genome
    phage_contig=$(awk -F'\t' -v phage="${phage}" '(($1 == phage) && ($3 == "yes")) {print $4}' ../induced_phages_2.txt)

    echo ${phage_contig} > phage_contig.id

    seqtk subseq \
        -l 80 \
        ../induced_phages.fasta \
        phage_contig.id > ${phage}_contig.fasta

    rm phage_contig.id



    ### trim the second copy of DTR (PhageTerm does not automatically trim it when doing sequence rearrangements)
    if [[ "${phage}" != "Tulp" ]]
    then
        seqtk trimfq -e 127 ${phage}_contig.fasta | seqtk seq -l 80 > tmp.fasta
        mv tmp.fasta ${phage}_contig.fasta
    fi



    ### reverse complement genomes so that most ORFs are on the forward strand
    if [[ ("${phage}" == "Pioen") || ("${phage}" == "Aster") || ("${phage}" == "Lelie") ]]
    then
        seqtk seq -l 80 -r ${phage}_contig.fasta > tmp.fasta
        mv tmp.fasta ${phage}_contig.fasta
    fi



    ### get prophage-containing contig
    host_strain=$(awk -F'\t' -v phage="${phage}" '(($1 == phage) && ($3 == "yes")) {print $8}' ../induced_phages_2.txt)

    if [[ "${phage}" == "Roos" ]]
    then

        host_contig='NODE_1_length_167466_cov_21.860236'

    else

        CP=$(awk -F'\t' -v phage="${phage}" '(($1 == phage) && ($3 == "yes")) {print $2}' ../induced_phages_2.txt)

        host_contig=$(awk -F'\t' -v CP="${CP}" '$1 == CP {print $3}' ../../CPs_refine_coo/prophages_refined_coo.txt)

    fi

    echo ${host_contig} > host_contig.id

    seqtk subseq \
        -l 80 \
        ../../Fprau_assemblies/${host_strain}_contigs.fasta \
        host_contig.id > ${phage}_host_contig.fasta

    rm host_contig.id

done



module purge; module load Anaconda3; module list
conda activate PhageTerm_env; conda list

for phage in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie'
do

    sample=$(awk -F'\t' -v phage="${phage}" '(($1 == phage) && ($3 == "yes")) {print $7}' ../induced_phages_2.txt)


    DD=../../Induction_raw_reads/X204SC23012117-Z01-F004/01.RawData/${sample}/
    F1=$(ls ${DD}/${sample}_*_1.fq.gz)
    F2=$(ls ${DD}/${sample}_*_2.fq.gz)

    cp -p ${F1} .
    cp -p ${F2} .

    B1=$(basename ${F1} .fq.gz)
    B2=$(basename ${F2} .fq.gz)

    gzip -d ${B1}.fq.gz
    gzip -d ${B2}.fq.gz


    python ${HOME}/SOFTWARE/phageterm-1.0.12/PhageTerm.py \
        --phagename ${phage} \
        --ref ${phage}_contig.fasta \
        --fastq  ${B1}.fq \
        --paired ${B2}.fq \
        --host ${phage}_host_contig.fasta \
        --core ${SLURM_CPUS_PER_TASK}


    rm ${B1}.fq ${B2}.fq ${phage}_host_contig.fasta
    chmod 440 ${phage}_*

done
