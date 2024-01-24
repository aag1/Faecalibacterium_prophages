sessionInfo()



STRAIN_ID <- commandArgs(trailingOnly = T)[1]


f1 <- paste0('../Fprau_map_reads/BED_FILES/', STRAIN_ID, '_contigs.bed')
t1 <- read.table(f1, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)


f2 <- paste0('../Fprau_map_reads/mapped_own_reads/', STRAIN_ID, '_own_reads.depth.txt')
t2 <- read.table(f2, sep = '\t', header = F, stringsAsFactors = F)


Wq <- list()
Q <- list()



for (CONTIG_ID in rownames(t1)[t1[, 2] > 50000]) {

    if (!(CONTIG_ID %in% t2[, 1])) { next }


    len <- t1[CONTIG_ID, 2]


    i <- which(t2[, 1] == CONTIG_ID)
    coo <- t2[i, 2]
    cov <- t2[i, 3]


    COVERAGE <- rep(0, len)
    COVERAGE[coo] <- cov


    w_center <- seq(from = 1501, to = len - 1500, by = 500)
    w_from <- w_center - 1500
    w_to <- w_center + 1500
    w_to[ length(w_to) ] <- len


    # sliding window centers
    Wq[[ CONTIG_ID ]] <- w_center


    # depth of coverage per window
    Q[[ CONTIG_ID ]] <- sapply(seq_along(w_center), function (i) {

                mean( COVERAGE[ w_from[i]:w_to[i] ] )

    })

}



save(Wq, Q, file = paste0(STRAIN_ID, '_own_reads_coverage_summary.RData'))
