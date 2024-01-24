library(bio3d)
source('function_plot_msa.R')
sessionInfo()



pdf('DGR_repeats_ali.pdf', height = 7, width = 7.5)

par(mar = c(1, 6, 2, 1), cex.axis = 0.55, mgp = c(3, 0.5, 0), tck = -0.01, cex = 0.75)

plot(
    NA,
    xlim = c(0, 137),
    ylim = c(0, 30),
    xaxs = 'i',
    yaxs = 'i',
    axes = F,
    ann = F
)



nt_col <- setNames(c('lightgoldenrod1', 'darkseagreen3', 'cornflowerblue', 'lightcoral', 'white'), c('A', 'T', 'G', 'C', '-'))

rect(xleft = 106, xright = 137, ybottom = 28, ytop = 30)

xL <- 109

for (n in names(nt_col)[1:4]) {

    rect(xleft = xL, xright = xL + 1, ybottom = 28.5, ytop = 29.5, col = nt_col[n], border = NA)

    text(n, x = xL + 2, y = 29, adj = 0)

    xL <- xL + 7

}



yBottom <- 0

for (x in rev(c('Tulp', 'Roos', 'Pioen', 'Lelie', 'CP2', 'CP11'))) {

    f <- paste0(x, '_DGR_repeats_ali.fasta')


    A <- read.fasta(f)$ali
    A <- toupper(A)


    plot_msa(A, nt_col, xLeft = 0, yBottom = yBottom)


    yBottom <- yBottom + nrow(A) + 3

}

dev.off()
