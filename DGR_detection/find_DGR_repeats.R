.libPaths(paste0(Sys.getenv('HOME'), '/SOFTWARE/R_LIB'))
library(IRanges)
sessionInfo()




args <- commandArgs(trailingOnly = T)

p <- args[1]


blastn_file <- paste0('BLASTN_out/', p, '_repeats.6a.txt')

t <- read.table(blastn_file, sep = '\t', header = F, stringsAsFactors = F)

colnames (t) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'qcovs', 'sstart', 'send', 'slen', 'sstrand', 'evalue', 'bitscore', 'nident', 'qseq', 'sseq')


proteome_file <- paste0(ifelse(grepl('^CP', p), '../CPs_proteomes/', '../Phages_proteomes/'), p, '_proteome.txt')

q <- read.table(proteome_file, sep = '\t', header = F, stringsAsFactors = F)

colnames (q) <- c('orf_name', 'orf_from', 'orf_to', 'orf_strand')




# ------------------------------ find repeats ------------------------------ #
bad <- c()

for (i in 1:nrow(t)) {

    ### exclude hits between a whole genome and itself
    if (t$qlen[i] == t$nident[i]) { bad <- c(bad, i) }


    ### exclude duplicates
    if (i > 1) {

        for (j in 1:(i-1)) {

            if (t$sstrand[i] != t$sstrand[j]) { next }

            if ((t$sstrand[i] == 'plus') &
                (t$qstart[i] == t$sstart[j]) & (t$qend[i] == t$send[j]) &
                (t$sstart[i] == t$qstart[j]) & (t$send[i] == t$qend[j])) { bad <- c(bad, i) }

            if ((t$sstrand[i] == 'minus') &
                (t$qstart[i] == t$send[j]) & (t$qend[i] == t$sstart[j]) &
                (t$sstart[i] == t$qend[j]) & (t$send[i] == t$qstart[j])) { bad <- c(bad, i) }

        }

    }


    ### exclude: short alignments, low identity alignments, (almost) identical alignments
    if ((t$length[i] < 90) | (t$pident[i] < 50) | (t$mismatch[i] < 5)) { bad <- c(bad, i) }

}

if (length(bad) > 0) { t <- t[-bad, ] }




# ------------------------------ find repeats with A mismatches ------------------------------ #
DF <- data.frame(NULL)

if (nrow(t) > 0) {

    for (i in 1:nrow(t)) {

        qseq <- strsplit(t$qseq[i], '')[[1]]

        sseq <- strsplit(t$sseq[i], '')[[1]]


        idx <- which(qseq != sseq)


        for (qN in c('A', 'T')) {

            ratio <- sum(qseq[idx] == qN) / length(idx)

            qTR <- (ratio >= 0.75)

            if (qTR) { break }

        }


        for (sN in c('A', 'T')) {

            ratio <- sum(sseq[idx] == sN) / length(idx)

            sTR <- (ratio >= 0.75)

            if (sTR) { break }

        }


        if ((!qTR) & (!sTR)) { next }

        if (qTR & sTR) { stop('A pair of repeats appears to have two TRs!') }

        if (qTR) {

            df <- data.frame(
                phage = p,
                TR_coo1 = t$qstart[i],
                TR_coo2 = t$qend[i],
                VR_coo1 = t$sstart[i],
                VR_coo2 = t$send[i],
                TR_nt = qN,
                stringsAsFactors = F
            )

        }

        if (sTR) {

            df <- data.frame(
                phage = p,
                TR_coo1 = t$sstart[i],
                TR_coo2 = t$send[i],
                VR_coo1 = t$qstart[i],
                VR_coo2 = t$qend[i],
                TR_nt = sN,
                stringsAsFactors = F
            )

        }


        DF <- rbind(DF, df)

    }

}




# ------------------------------ ORF overlaps ------------------------------ #
if (nrow(DF) > 0) {

    DF <- DF[order(DF$VR_coo1, decreasing = F), ]


    DF[, c('TR_orf', 'TR_pos', 'VR_orf', 'VR_pos', 'VR_aa1', 'VR_aa2', 'VR_aaL')] <- NA


    ir1 <- IRanges(start = q$orf_from, end = q$orf_to)


    for (i in 1:nrow(DF)) {

        for (x in c('TR', 'VR')) {

            coo1 <- DF[i, paste0(x, '_coo1')]
            coo2 <- DF[i, paste0(x, '_coo2')]
            ir2 <- IRanges(start = coo1, end = coo2)


            ### overalps
            OV <- findOverlaps(ir2, ir1)
            OV <- as.data.frame(OV)

            idx <- OV$subjectHits
            if (length(idx) == 0) { next }

            DF[i, paste0(x, '_orf')] <- paste(q$orf_name[idx], collapse = ';')


            ### overlap type
            OV_type <- sapply(idx, function (j) {

                boo <- ((q$orf_from[j] <= coo1) & (coo2 <= q$orf_to[j]))

                ifelse(boo, 'belongs', 'overlaps')

            })

            DF[i, paste0(x, '_pos')] <- paste(OV_type, collapse = ';')


            ### VR target protein info
            if (x == 'VR') {

                if (length(idx) > 1) { stop('VR overlaps with multiple ORFs!') }


                if (q$orf_strand[idx] == '+') {

                    DF[i, 'VR_aa1'] <- ceiling((coo1 - q$orf_from[idx] + 1) / 3)
                    DF[i, 'VR_aa2'] <- ceiling((coo2 - q$orf_from[idx] + 1) / 3)

                } else {

                    DF[i, 'VR_aa1'] <- ceiling((q$orf_to[idx] - 3 - coo2 + 1) / 3)
                    DF[i, 'VR_aa2'] <- ceiling((q$orf_to[idx] - 3 - coo1 + 1) / 3)

                }


                DF[i, 'VR_aaL'] <- (q$orf_to[idx] - q$orf_from[idx] + 1) / 3 - 1

            }

        }

    }


    write.table(DF, sep = '\t', quote = F, row.names = F, file = paste0(p, '_DGR_repeats.txt'))

}
