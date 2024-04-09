sessionInfo()



t <- read.table('PhIFM7_102A3_split_coo.txt', sep = '\t', header = T, stringsAsFactors = F)

b <- read.table('../Fprau_map_reads/BED_FILES/FM7_S1_contigs.bed', sep = '\t', row.names = 1, header = F, stringsAsFactors = F)
b <- b[b[, 2] > 50000, ]



pdf('Mushu_transpos_sites.pdf', height = 5, width = 7.5)

layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = T), height = c(3.5, 2), width = c(3.25, 2))

par(las = 1, mgp = c(3, 0.75, 0))

COL <- adjustcolor('blueviolet', alpha.f = 0.5)




# ------------------------------ panel A ------------------------------ #
par(mar = c(5, 5, 1.25, 1))

d <- table(t$pos2_adj[t$cigar2_type == 'SM'])

k <- barplot(d, xlab = '', ylab = 'Number of VLP reads', names.arg = '', col = COL)

text(names(d), x = k, y = -0.025 * par()$usr[4], srt = 90, adj = 1, cex = 0.5, xpd = T)

mtext('NODE_11 nucleotide fused at 5\'-end', side = 1, line = 2, cex = 4/5)

mtext('A', side = 2, line = 3.5, at = par()$usr[4], cex = 1.75)




# ------------------------------ panel B ------------------------------ #
par(mar = c(5, 5, 1.25, 1))

d <- table(t$pos2_adj[t$cigar2_type == 'MS'])

k <- barplot(d, xlab = '', ylab = 'Number of VLP reads', names.arg = '', col = COL)

text(names(d), x = k, y = -0.025 * par()$usr[4], srt = 90, adj = 1, cex = 0.5, xpd = T)

mtext('NODE_11 nucleotide fused at 3\'-end', side = 1, line = 2, cex = 4/5, adj = 0.9)

mtext('B', side = 2, line = 3.5, at = par()$usr[4], cex = 1.75)




# ------------------------------ panel C ------------------------------ #
t <- t[t$pos2_adj == 53400, ]


par(mar = c(2, 5, 1, 1))

plot(
    NA,
    xlim = c(0, sum(b[, 2])),
    ylim = c(0, 12),
    xlab = '',
    ylab = 'Mushu transposition sites',
    xaxs = 'i',
    yaxs = 'i',
    bty = 'n',
    axes = F
)

axis(side = 2, at = c(0, 5, 10))

mtext('contig', side = 1, line = 0.5, at = -25000, adj = 1, cex = 0.75)


Xleft <- 0

for (CONTIG_ID in rownames(b)) {

    len <- b[CONTIG_ID, 2]

    w_center <- seq(from = 1501, to = len - 1500, by = 500)
    w_from <- w_center - 1500
    w_to <- w_center + 1500
    w_to[ length(w_to) ] <- len

    V <- unique(t$pos1_adj[t$contig1 == CONTIG_ID])
    Q <- sapply(seq_along(w_center), function (i) { sum((V >= w_from[i]) & (V <= w_to[i])) })

    lines(Xleft + w_center, Q, col = COL, xpd = T)

    rect(xleft = Xleft, xright = Xleft + len, ybottom = par()$usr[3], ytop = par()$usr[4])

    n <- sub('^NODE_([0-9]+)_.+$', '\\1', CONTIG_ID)
    mtext(n, side = 1, line = 0.5, at = Xleft + len / 2, cex = 0.75, adj = ifelse(n == '13', 0, 0.5))

    Xleft <- Xleft + len

}


mtext('C', side = 2, line = 3.5, at = par()$usr[4] * 1.2, cex = 1.75)


dev.off()
