#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



wget --no-verbose ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/Pfam-A.hmm.gz

gzip -d Pfam-A.hmm.gz



echo -e 'NAME\tACC\tDESC' > PfamA_36.0_domains.txt
grep -A2 '^NAME ' Pfam-A.hmm |
sed -E 's/^[A-Z]+ +//' |
tr '\n' '\t' |
sed 's/\t--\t/\n/g' |
sed '$ s/\t$/\n/' >> PfamA_36.0_domains.txt



chmod 440 Pfam-A.hmm PfamA_36.0_domains.txt
