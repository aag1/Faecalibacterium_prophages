sessionInfo()




pdf('phage_ecology.pdf', height = 2.25, width = 7.5)

layout(matrix(1:3, nrow = 1, byrow = T), width = c(3.6, 4.3, 2.1))




# ------------------------------ panel A ------------------------------ #
phages <- c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie')
hosts <- c(paste0('s__Faecalibacterium prausnitzii', c('', '_C', '_F', '_G', '_H')), 's__Faecalibacterium sp900539885')
M <- matrix(0, nrow = 6, ncol = 5, dimnames = list(hosts, phages))

Discovery <- M
Discovery['s__Faecalibacterium prausnitzii_C', 'Tulp']  <- 1
Discovery['s__Faecalibacterium prausnitzii_C', 'Roos']  <- 1
Discovery['s__Faecalibacterium prausnitzii', 'Pioen']   <- 1
Discovery['s__Faecalibacterium prausnitzii_F', 'Pioen'] <- 1
Discovery['s__Faecalibacterium sp900539885', 'Pioen']   <- 1
Discovery['s__Faecalibacterium prausnitzii_F', 'Aster'] <- 1
Discovery['s__Faecalibacterium prausnitzii_G', 'Lelie'] <- 1

Infection <- M
t <- read.table('../Phages_hosts/blast_hosts_summary.txt', sep = '\t', header = T, stringsAsFactors = F)
for (i in 1:nrow(t)) {
    Infection[t$host_taxo[i], t$phage[i]] <- 1
}

CRISPR <- M
t <- read.table('../Phages_hosts/crispr_hosts_summary.txt', sep = '\t', header = T, stringsAsFactors = F)
for (i in 1:nrow(t)) {
    CRISPR[t$host_taxo[i], t$phage[i]] <- 1
}


par(mar = c(4, 8, 4, 1), las = 1)


plot(NA, xlim = c(0, length(phages)), ylim = c(length(hosts), 0), xaxs = 'i', yaxs = 'i', ann = F, axes = F)

for (i in 1:nrow(M)) {

    for (j in 1:ncol(M)) {

        if (Discovery[i, j] == 1) { points(x = j - 0.5, y = i - 0.5, pch = 15, col = 'grey75', cex = 3) }

        if (Infection[i, j] == 1) { points(x = j - 0.5, y = i - 0.5, pch = 4, cex = 1.5) }

        if (CRISPR[i, j] == 1) { points(x = j - 0.5, y = i - 0.5, pch = 16, cex = 1.5) }

    }

}

hosts <- sub('^s__Faecalibacterium', 'F.', hosts)
hosts <- gsub('_', ' ', hosts)
text(labels = hosts, x = -0.5, y = seq_along(hosts) - 0.5, adj = 1, xpd = T)

phages[phages == 'Tulp'] <- 'Mushu'
text(labels = phages, x = seq_along(phages) - 0.5, y = -0.6, srt = 90, adj = 0, xpd = T)


points(x = -2.75, y = 7.5, pch = 15, col = 'grey75', cex = 3, xpd = T)
text('Detection', x = -2.25, y = 7.5, cex = 0.8, adj = 0, xpd = T)
points(x = -0.25, y = 7.5, pch = 4, cex = 1.5, xpd = T)
text('Detection UHGG', x = 0, y = 7.5, cex = 0.8, adj = 0, xpd = T)
points(x = 2.75, y = 7.5, pch = 16, cex = 1.5, xpd = T)
text('CRISPR UHGG', x = 3, y = 7.5, cex = 0.8, adj = 0, xpd = T)
rect(xleft = -3.25, xright = 5.25, ybottom = 8.25, ytop = 6.75, xpd = T)




# ------------------------------ panel B ------------------------------ #
DF <- read.table('LLD_IBD_phages_logistic_regression.txt', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

DF$stars <- sapply(DF$P_adj, function (x) {

    stars <- ''
    if (x < 0.05) { stars <- '*' }
    if (x < 0.01) { stars <- '**' }
    if (x < 0.001) { stars <- '***' }
    return(stars)

})

DF$LLD_pos <- DF$LLD_pos * 100
DF$IBD_pos <- DF$IBD_pos * 100

M <- as.matrix(t(DF[c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie'), c('LLD_pos', 'IBD_pos')]))
colnames(M)[colnames(M) == 'Tulp'] <- 'Mushu'


par(mar = c(3, 8, 2, 1), las = 1)


bp <- barplot(
    height = M,
    beside = T,
    ylim = c(0, 40),
    ylab = 'Positive samples, %'
)

for (i in c(1:2, 4:5)) { lines(x = bp[, i], y = rep(max(M[, i]) + 2, 2), xpd = T) }

text(DF$stars[1:5], x = apply(bp, 2, sum) / 2, y = apply(M, 2, max) + 4, cex = 1.5, xpd = T)


legend(
    'topright',
    c('LLD', '1000IBD'),
    fill = gray.colors(2),
    ncol = 1,
    inset = c(0, -0.1),
    xpd = T
)


mtext(c('A', 'B'), side = 2, line = c(26, 4), at = par()$usr[4], cex = 2)




# ------------------------------ panel C ------------------------------ #
M <- as.matrix(t(DF['Host', c('LLD_pos', 'IBD_pos')]))


bp <- barplot(
    height = M,
    beside = T,
    xlim = c(0.5, 3.1),
    ylim = c(0, 100),
    names.arg = '',
    ylab = 'Positive samples, %'
)

lines(x = bp, y = rep(max(M) + 5, 2), xpd = T)

text(DF$stars[6], x = apply(bp, 2, sum) / 2, y = apply(M, 2, max) + 10, cex = 1.5, xpd = T)

mtext(
    'Faecalibacterium',
    side = 1,
    line = 1,
    at = 0,
    cex = 2/3,
    font = 3
)


mtext('C', side = 2, line = 4, at = par()$usr[4], cex = 2)


dev.off()
