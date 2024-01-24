source('plot_panel_functions.R')
sessionInfo()




# ------------------------------ read data ------------------------------ #
STRAIN_ID <- 'HTF-128_S3'


file1 <- paste0('../Fprau_proteomes/', STRAIN_ID, '_proteome.gff')
P <- read.table(file1, sep = '\t', header = F, stringsAsFactors = F)[, c(1, 4, 5, 7)]


X <- list()
for (n in c('ct2', 'vs2', 'geNomad')) {

    f <- paste0('../Fprau_prophage_tools/Fprau_prophages_', n, '.txt')
    t <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)
    X[[ n ]] <- t[t$strain_id == STRAIN_ID, ]

}


file2 <- paste0('../Fprau_map_reads/BED_FILES/', STRAIN_ID, '_contigs.bed')
tab <- read.table(file2, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)
tab <- tab[tab[, 2] > 50000, ]


load(paste0('CONTIGS_INFO/', STRAIN_ID, '_vlp_coverage_summary.RData'))          # W, B, D

load(paste0('CONTIGS_INFO/', STRAIN_ID, '_own_reads_coverage_summary.RData'))    # Wq, Q

load(paste0('CONTIGS_INFO/', STRAIN_ID, '_nt_content.RData'))                    # Wn, N


DF <- read.table('prophages_raw_coo.txt', sep = '\t', header = T, stringsAsFactors = F)

study_col <- setNames(c('lightpink2', 'steelblue2'), c('VanEspen_2021', 'Shkoporov_2019'))

nt_col <- setNames(c('lightgoldenrod1', 'darkseagreen', 'cornflowerblue', 'lightcoral'), c('A', 'T', 'G', 'C'))




# ------------------------------ initiate plot ------------------------------ #
pdf('CP11_example.pdf', height = 5.5, width = 7.5)

layout(matrix(1:6, nrow = 6), height = c(0.9, 1, 0.8, 1, 1, 0.8))

par(oma = c(3, 1, 1, 0), las = 1)

Xmax <- 216000

Q_range <- c(40, 100)

CONTIG_ID <- DF$contig_id[DF$strain_id == STRAIN_ID]

CONTIG_ID_short <- sub('^(NODE_[0-9]+)_.+$', '\\1', CONTIG_ID)

len <- tab[CONTIG_ID, 2]
LEN <- ceiling(len / 1000) * 1000

Xaxis <- seq(from = 0, to = LEN, by = 50000)
if (LEN - Xaxis[length(Xaxis)] < 25000) { Xaxis[length(Xaxis)] <- LEN } else { Xaxis <- c(Xaxis, LEN) }

df <- DF[(DF$strain_id == STRAIN_ID) & (DF$contig_id == CONTIG_ID), ]




# ------------------------------ plot ORFs ------------------------------ #
par(mar = c(1.5, 6, 1.5, 1))

plot_ORFs_panel(Xmax, df, CONTIG_ID, len, cex_lab = 0.9)

mtext('A', side = 2, line = 5, at = par()$usr[4], cex = 1.5)




# ------------------------------ plot nt content ------------------------------ #
plot_nt_content_panel(Xmax, Xaxis, df, CONTIG_ID, Wn, N, nt_col, x_axis_labs = F)

legend(
    lty = 1,
    col = nt_col,
    legend = names(nt_col),
    'topleft', inset = c(0.01, -0.3), xpd = T, ncol = 4,
    cex = 0.75
)

mtext('B', side = 2, line = 5, at = par()$usr[4], cex = 1.5)




# ------------------------------ plot prophages predicted by tools ------------------------------ #
plot_tools_panel(Xmax, Xaxis, df, CONTIG_ID, X, len, cex_lab = 0.75, lwd_signal = 6, x_axis_labs = F)

mtext('C', side = 2, line = 5, at = par()$usr[4], cex = 1.5)




# ------------------------------ plot vlp coverage breadth ------------------------------ #
vlp_cov_breadth_panel(Xmax, Xaxis, df, CONTIG_ID, W, B, study_col, x_axis_labs = F)

mtext('D', side = 2, line = 5, at = par()$usr[4], cex = 1.5)




# ------------------------------ plot vlp coverage depth ------------------------------ #
vlp_cov_depth_panel(Xmax, Xaxis, df, CONTIG_ID, W, D, study_col, 2, x_axis_labs = F)

legend(
    lty = 1,
    col = sapply(names(study_col), function (x) study_col[ x ]),
    legend = sub('VanEspen', 'Van Espen', sub('_', ' et al. ', names(study_col))),
    title = 'VLP reads:',
    'topleft', inset = c(0.01, -0.3), xpd = T, ncol = 1,
    cex = 0.75
)

mtext('E', side = 2, line = 5, at = 100, cex = 1.5)




# ------------------------------ plot own reads coverage depth ------------------------------ #
own_reads_panel(Xmax, Xaxis, df, CONTIG_ID, Wq, Q, Q_range)

legend(
    lty = 1,
    col = 'gray25',
    legend = 'own reads',
    'topleft', inset = c(0.01, -0.4), xpd = T,
    cex = 0.75
)

mtext(paste0(sub('_S[0-9]+$', '', STRAIN_ID), '  ', CONTIG_ID_short, ', kb'), side = 1, line = 2.75, at = max(Xaxis) / 2, cex = 0.9)

mtext('F', side = 2, line = 5, at = par()$usr[4], cex = 1.5)

dev.off()
