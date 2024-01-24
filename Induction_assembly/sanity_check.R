library(seqinr)



S <- read.fasta('induced_phages.fasta')

tab <- read.table('induced_phages_2.txt', sep = '\t', header = T, stringsAsFactors = F)

coo <- read.table('integr_sites_coo.txt', sep = '\t', header = T, stringsAsFactors = F)



### Were the pre-PhageTerm rearrangements (second DTR copy trimmed, reverse complement) made as expected?
for (x in c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie')) {

    contig <- tab$contig_name[ (tab$phage_name == x) & (tab$repres == 'yes') ]

    s0 <- S[[ contig ]]

    s1 <- read.fasta(paste0('PhageTerm_OUT/', x, '_contig.fasta'))[[1]]


    if (x != 'Tulp') { s0 <- s0[1:(length(s0) - 127)] }

    if (x %in% c('Pioen', 'Aster', 'Lelie')) { s0 <- rev(comp(s0, ambiguous = T)) }


    print(all(s0 == s1))

}



### Were the PhageTerm rearrangements made as expected?
DF <- data.frame(phage = c('Roos', 'Pioen', 'Aster', 'Lelie'), start = c(28477, 20612, 17572, 44778), stringsAsFactors = F)
DF$len <- NA

for (i in 1:nrow(DF)) {

    s1 <- read.fasta(paste0('PhageTerm_OUT/', DF$phage[i], '_contig.fasta'))[[1]]

    s2 <- read.fasta(paste0('PhageTerm_OUT/', DF$phage[i], '_sequence.fasta'))[[1]]


    L <- length(s1); DF$len[i] <- L

    N <- DF$start[i]

    x <- s1[c(N:L, 1:(N-1))]


    print(all(x == s2))

}



### Is the fusion made by PhageTerm present in the corresponding prophage?
DF$fusion <- DF$len - DF$start + 1

for (i in 1:nrow(DF)) {

    j <- which(coo$phage_name == DF$phage[i])[1]


    prophage_seq <- read.fasta(paste0('../Fprau_assemblies/', coo$strain[j], '_contigs.fasta'))[[ coo$contig[j] ]][ coo$coo1[j]:coo$coo2[j] ]

    if (i %in% c(1, 4)) { prophage_seq <- rev(comp(prophage_seq, ambiguous = T)) }

    prophage_seq <- paste(prophage_seq, collapse = '')


    s2 <- read.fasta(paste0('PhageTerm_OUT/', DF$phage[i], '_sequence.fasta'))[[1]]

    fusion_seq <- paste(s2[DF$fusion[i] + (-49):50], collapse = '')


    print(grepl(fusion_seq, prophage_seq))

}
