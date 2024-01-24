library(seqinr)
sessionInfo()



args <- commandArgs(trailingOnly = T)
seq1_fasta <- args[1]
seq2_fasta <- args[2]
blast_res <- args[3]
out_fasta <- args[4]



# ------------------------------ read data ------------------------------ #
S1 <- read.fasta(
    seq1_fasta,
    seqtype = 'DNA',
    forceDNAtolower = F
)


S2 <- read.fasta(
    seq2_fasta,
    seqtype = 'DNA',
    forceDNAtolower = F
)


tab <- read.table(
    blast_res,
    sep = '\t',
    header = F,
    stringsAsFactors = F
)

colnames(tab) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'qcovs', 'sstart', 'send', 'slen', 'sstrand', 'evalue', 'bitscore', 'nident')

tab <- tab[tab$length > 400, ]



# ------------------------------ process data ------------------------------ #
START <- tab$sstart[1]

S2[[1]] <- S2[[1]][ c(START:length(S2[[1]]), 1:(START-1)) ]

L <- c(S1, S2)

write.fasta(
    sequences = L,
    names = names(L),
    nbchar = 80,
    file.out = out_fasta
)
