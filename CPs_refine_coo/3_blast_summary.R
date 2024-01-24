.libPaths(paste0(Sys.getenv('HOME'), '/SOFTWARE/R_LIB'))
library(IRanges)
sessionInfo()



DATA <- data.frame(NULL)

for (n in 1:15) {

    for (w in c('viral_refseq', 'IMGVR')) {

        t <- read.table(paste0('CP', n, '_raw_vs_', w, '.txt'), sep = '\t', header = F, stringsAsFactors = F)

        colnames(t) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'qcovs', 'sstart', 'send', 'slen', 'sstrand', 'evalue', 'bitscore', 'nident')


        DF <- data.frame(NULL)


        for (x in unique(t$sseqid)) {

            idx <- which(t$sseqid == x)


            d <- t[idx, ]

            ir <- IRanges(
                start = sapply(1:nrow(d), function (i) ifelse(d$sstrand[i] == 'plus', d$sstart[i], d$send[i])),
                end   = sapply(1:nrow(d), function (i) ifelse(d$sstrand[i] == 'plus', d$send[i], d$sstart[i]))
            )

            cl <- reduce(ir)

            cl <- as.data.frame(cl)

            scovs <- round(sum(cl$width) / d$slen[1] * 100)


            df <- data.frame(
                prophage = paste0('CP', n),
                db_name = w,
                db_seqid = x,
                qcovs = d$qcovs[1],
                qcoo = paste0(min(d$qstart), '-', max(d$qend)),
                qlen = d$qlen[1],
                scovs,
                slen = d$slen[1],
                stringsAsFactors = F
            )


            DF <- rbind(DF, df)

        }


        DF <- DF[order(DF$scovs, decreasing = T), ]


        DATA <- rbind(DATA, DF)

    }

}

DATA <- DATA[DATA$scovs > 25, ]



write.table(DATA, sep = '\t', row.names = F, quote = F, file = 'blast_summary.txt')
