source('../CPs_refine_coo/function_dp.R')
sessionInfo()




# ------------------------------ read data ------------------------------ #
PHAGES <- c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie')

CPS <- paste0('CP', c(2, 4, 5, 6, 7, 11))



tab <- read.table('pairwise_identity.txt', sep = '\t', header = F, stringsAsFactors = F)

colnames(tab) <- c('seq1', 'seq2', 'match', 'length', 'ident')

p <- read.table('../CPs_refine_coo/prophages_refined_coo.txt', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

r <- read.table('../Viral_RefSeq_220/viral_refseq_taxo.txt', sep = '\t', header = T, row.names = 1, quote = '', fill = T, stringsAsFactors = F)

get_len <- function (x) {

    if (x %in% PHAGES) {

        l <- read.table(paste0('../Phages_map_reads/BED_FILES/', x , '_genome.bed'), sep = '\t', header = F, stringsAsFactors = F)[1, 3]

    } else if (x %in% CPS) {

        l <- p[x, 'new_to'] - p[x, 'new_from'] + 1

    } else {

        l <- r[x, 'genome_length']

    }

    return(l)

}

tab$seq1_len <- sapply(tab$seq1, get_len)
tab$seq2_len <- sapply(tab$seq2, get_len)



info <- data.frame(
    seq_id = c('NC_047911.1', 'NC_047913.1', 'NC_047914.1', 'NC_047915.1'),
    seq_name = c('Lagaffe*', 'Mushu*', 'Taranis', 'Toutatis'),
    stringsAsFactors = F
)

row.names(info) <- info$seq_id




# ------------------------------ plot ------------------------------ #
pdf('dotplots_VCs.pdf', height = 2, width = 7.5)

layout(matrix(1:4, nrow = 1, byrow = T))

par(mar = c(5, 4, 1, 1))

Xmax <- max(tab$seq2_len)
Ymax <- max(tab$seq1_len)

for (i in 1:nrow(tab)) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(0, Ymax),
        xaxs = 'i',
        asp = 1,
        axes = F,
        ann = F
    )


    dp(
        file = paste0('dotplots_data/', tab$seq1[i], '___', tab$seq2[i], '.gff'),
        yID = tab$seq1[i]
    )


    Xaxis <- seq(from = 0, floor(tab$seq2_len[i] / 10000) * 10000, by = 10000)
    axis(side = 1, pos = 0, at = Xaxis, labels = Xaxis / 1000, cex.axis = 0.75, tck = -0.025, mgp = c(3, 0.1, 0))


    Xlab <- info$seq_name[ info$seq_id == tab$seq2[i] ]
    mtext(Xlab, side = 1, line = 1.25, at = tab$seq2_len[i] / 2)

    pct <- paste0(format(tab$ident[i], nsmall = 1), '%')
    mtext(pct, side = 1, line = 2.75, at = tab$seq2_len[i] / 2, cex = 0.75)


    Yaxis <- seq(from = 0, floor(tab$seq1_len[i] / 10000) * 10000, by = 10000)
    axis(side = 2, pos = 0, at = Yaxis, labels = Yaxis / 1000, las = 1, cex.axis = 0.75, tck = -0.025, mgp = c(3, 0.5, 0))


    mtext(tab$seq1[i], side = 2, line = 1.75, at = tab$seq1_len[i] / 2)


    rect(xleft = 0, xright = tab$seq2_len[i], ybottom = 0, ytop = tab$seq1_len[i], xpd = T)

}

dev.off()
