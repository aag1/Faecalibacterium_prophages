library(seqinr)
sessionInfo()



t <- read.table('integr_sites_coo.txt', sep = '\t', header = T, stringsAsFactors = F)



DF <- data.frame(NULL)

for (i in 1:nrow(t)) {

    phage <- t$phage_name[i]
    if (phage == 'Roos') { next }

    S <- read.fasta(paste0('../Fprau_assemblies/', t$strain[i], '_contigs.fasta'), seqtype = 'DNA', forceDNAtolower = F)

    flank1 <- paste(S[[ t$contig[i] ]][ t$coo1[i] - 50:1 ], collapse = '')
    flank2 <- paste(S[[ t$contig[i] ]][ t$coo2[i] + 1:50 ], collapse = '')

    site <- paste(flank1, '|', flank2)

    df <- data.frame(phage, site, stringsAsFactors = F)
    DF <- rbind(DF, df)

}



i1 <- which(t$phage_name == 'Roos')[1]
i2 <- which(t$phage_name == 'Roos')[2]

S <- read.fasta('../Fprau_assemblies/FM8_S2_contigs.fasta', seqtype = 'DNA', forceDNAtolower = F)

flank1 <- paste(S[[ t$contig[i2] ]][ t$coo1[i2] - 50:1 ], collapse = '')
flank2 <- paste(S[[ t$contig[i1] ]][ t$coo2[i1] + 1:50 ], collapse = '')

site <- paste(flank1, '|', flank2)

df <- data.frame(phage = 'Roos', site, stringsAsFactors = F)
DF <- rbind(DF, df)



write.table(DF, sep = '\t', row.names = F, quote = F, file = 'integr_sites_seq.txt')
