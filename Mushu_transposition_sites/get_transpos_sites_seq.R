library(seqinr)
sessionInfo()



t <- read.table('PhIFM7_102A3_split_coo.txt', sep = '\t', header = T, stringsAsFactors = F)
t <- t[t$pos2_adj == 53400, ]
t <- t[, c('contig1', 'pos1_adj')]
t$n_reads <- sapply(1:nrow(t), function (i) sum((t$contig1 == t$contig1[i]) & (t$pos1_adj == t$pos1_adj[i])))
t <- unique(t)

b <- read.table('../Fprau_map_reads/BED_FILES/FM7_S1_contigs.bed', sep = '\t', row.names = 1, header = F, stringsAsFactors = F)

S <- read.fasta('../Fprau_assemblies/FM7_S1_contigs.fasta', seqtype = 'DNA', forceDNAtolower = F)



L <- list()

for (i in 1:nrow(t)) {

    contig <- t$contig1[i]

    N <- t$pos1_adj[i]

    from <- max(1, N - 49)

    to <- min(N + 50, b[contig, 2])



    flank1 <- paste(S[[ contig ]][ from:N ], collapse = '')

    flank2 <- paste(S[[ contig ]][ (N+1):to ], collapse = '')

    t$site[i] <- paste(flank1, '|', flank2)



    L[[paste0(contig, '_', from, '_', to)]] <- S[[ contig ]][ from:to ]

}



write.table(t, sep = '\t', row.names = F, quote = F, file = 'Tulp_transpos_sites_seq.txt')
write.fasta(names = names(L), sequences = L, nbchar = 80, file.out = 'Tulp_transpos_sites_seq.fasta')
