library(bio3d)



args <- commandArgs(trailingOnly = T)
seq1 <- args[1]
seq2 <- args[2]



ali_fasta <- paste0(seq1, '_', seq2, '_ali.fasta')
ALI <- read.fasta(ali_fasta)$ali



MATCH <- sum(apply(ALI, 2, function (v) (v[1] == v[2]) & (v[1] %in% c('a', 't', 'g', 'c'))))

LENGTH <- ncol(ALI)



X <- MATCH / LENGTH * 100

X <- round(X, digits = 1)

cat(paste0(seq1, '\t', seq2, '\t', MATCH, '\t', LENGTH, '\t', X, '\n'))
