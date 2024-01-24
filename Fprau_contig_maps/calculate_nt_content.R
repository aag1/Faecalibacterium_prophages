library(seqinr)
sessionInfo()



STRAIN_ID <- commandArgs(trailingOnly = T)[1]


f1 <- paste0('../Fprau_map_reads/BED_FILES/', STRAIN_ID, '_contigs.bed')
t1 <- read.table(f1, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)


f2 <- paste0('../Fprau_assemblies/', STRAIN_ID, '_contigs.fasta')
S <- read.fasta(f2, forceDNAtolower = F)


Wn <- list()
N <- list()



for (CONTIG_ID in rownames(t1)[t1[, 2] > 50000]) {

    len <- t1[CONTIG_ID, 2]


    w_center <- seq(from = 1501, to = len - 1500, by = 500)
    w_from <- w_center - 1500
    w_to <- w_center + 1500
    w_to[ length(w_to) ] <- len
    w_size <- w_to - w_from + 1


    # sliding window centers
    Wn[[ CONTIG_ID ]] <- w_center


    # base content per window
    N[[ CONTIG_ID ]] <- lapply(c('A', 'T', 'G', 'C'), function (b) {

            V <- ifelse(S[[ CONTIG_ID ]] == b, 1, 0)

            sapply(seq_along(w_center), function (i) {

                sum( V[ w_from[i]:w_to[i] ] ) / w_size[i] * 100

            })

    })

    names(N[[ CONTIG_ID ]]) <- c('A', 'T', 'G', 'C')

}



save(Wn, N, file = paste0(STRAIN_ID, '_nt_content.RData'))
