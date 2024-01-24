source('function_dp.R')
sessionInfo()



tab <- read.table('pairs_for_dotplots.txt', sep = '\t', header = T, stringsAsFactors = F)
tab$seq2 <- sub('^([^\\|]+)\\|.+$', '\\1', tab$seq2)
t <- read.table('new_coo.txt', sep = '\t', header = T, stringsAsFactors = F)



pdf('dotplots.pdf', height = 5.5, width = 7.5)
par(mar = c(3, 0, 0, 0), oma = c(0, 1, 0, 0))
layout(matrix(1:nrow(tab), ncol = 5, byrow = T))


Xmax <- max(tab$seq2_len)
Ymax <- max(tab$seq1_len)


for (i in 1:nrow(tab)) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(0, Ymax),
        yaxs = 'i',
        asp = 1,
        axes = F,
        ann = F
    )


    dp(
        file = paste0('dotplots_data/', tab$seq1[i], '___', tab$seq2[i], '.gff'),
        yID = tab$seq1[i]
    )


    Xaxis <- seq(from = 0, floor(tab$seq2_len[i] / 10000) * 10000, by = 10000)
    axis(side = 1, pos = 0, at = Xaxis, labels = Xaxis / 1000, cex.axis = 0.5, tck = -0.025, mgp = c(3, 0.1, 0))


    Xlab <- tab$seq2[i]
    if (tab$seq2_strand[i] == 'r') { Xlab <- paste(Xlab, 'r.c.') }
    mtext(Xlab, side = 1, line = 1, at = tab$seq2_len[i] / 2, cex = ifelse(tab$seq2_db[i] == 'viral_refseq', 0.5, 0.4))


    Yaxis <- seq(from = 0, floor(tab$seq1_len[i] / 10000) * 10000, by = 10000)
    axis(side = 2, pos = 0, at = Yaxis, labels = Yaxis / 1000, las = 1, cex.axis = 0.5, tck = -0.025, mgp = c(3, 0.5, 0))


    mtext(tab$seq1[i], side = 2, line = -0.25, at = tab$seq1_len[i] / 2, cex = 0.5)


    rect(xleft = 0, xright = tab$seq2_len[i], ybottom = 0, ytop = tab$seq1_len[i], xpd = T)


    j <- which(t$prophage == tab$seq1[i])
    lines(x = rep(tab$seq2_len[i] + 1500, 2), y = c(t$from[j], t$to[j]), lwd = 2, lend = 'butt', col = 'deepskyblue')
 
}


dev.off()
