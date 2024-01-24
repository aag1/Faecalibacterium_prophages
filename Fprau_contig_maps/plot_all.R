source('plot_panel_functions.R')
sessionInfo()




STRAIN_ID <- commandArgs(trailingOnly = T)[1]




# ------------------------------ read data ------------------------------ #
DF <- data.frame(NULL)
file0 <- 'prophages_raw_coo.txt'
if (file.exists(file0)) {

    DF <- read.table(file0, sep = '\t', header = T, stringsAsFactors = F)
    DF <- DF[DF$strain_id == STRAIN_ID, ]

}


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
if (nrow(tab) == 0) { stop('Strain ', STRAIN_ID, ' has no contigs > 50 kb!') }


load(paste0('CONTIGS_INFO/', STRAIN_ID, '_vlp_coverage_summary.RData'))          # W, B, D

load(paste0('CONTIGS_INFO/', STRAIN_ID, '_own_reads_coverage_summary.RData'))    # Wq, Q

load(paste0('CONTIGS_INFO/', STRAIN_ID, '_nt_content.RData'))                    # Wn, N


study_col <- setNames(c('lightpink2', 'steelblue2'), c('VanEspen_2021', 'Shkoporov_2019'))

nt_col <- setNames(c('lightgoldenrod1', 'darkseagreen', 'cornflowerblue', 'lightcoral'), c('A', 'T', 'G', 'C'))




# ------------------------------ initiate plot ------------------------------ #
pdf(paste0(STRAIN_ID, '_contig_maps.pdf'), height = 8.3, width = 11.7)

layout(matrix(1:6, nrow = 6), height = c(0.5, 1, 0.75, 1, 1, 1))

par(
    oma = par()$oma + c(2, 0, 0.5, 0),
    mgp = par()$mgp + c(1.25, 0, 0),
    las = 1,
    cex.lab = 1.5
)

Xmax <- 650000

Q_range <- range(unlist(lapply(Q, range)))
Q_range <- c(floor(Q_range[1] / 50) * 50, ceiling(Q_range[2] / 50) * 50)




for (CONTIG_ID in rownames(tab)) {

    CONTIG_ID_short <- sub('^(NODE_[0-9]+)_.+$', '\\1', CONTIG_ID)

    len <- tab[CONTIG_ID, 2]
    LEN <- ceiling(len / 1000) * 1000

    Xaxis <- seq(from = 0, to = LEN, by = 50000)
    if (LEN - Xaxis[length(Xaxis)] < 25000) { Xaxis[length(Xaxis)] <- LEN } else { Xaxis <- c(Xaxis, LEN) }

    df <- DF[DF$contig_id == CONTIG_ID, ]




    # ------------------------------ plot ORFs ------------------------------ #
    par(mar = c(0.2, 7, 2.1, 2.1))

    plot_ORFs_panel(Xmax, df, CONTIG_ID, len)




    # ------------------------------ plot nt content ------------------------------ #
    par(mar = c(3.1, 7, 2.1, 2.1))

    plot_nt_content_panel(Xmax, Xaxis, df, CONTIG_ID, Wn, N, nt_col)

    legend(
        lty = 1,
        col = nt_col,
        legend = names(nt_col),
        'topright', inset = c(-0.01, -0.15), xpd = T, ncol = 4
    )




    # ------------------------------ plot prophages predicted by tools ------------------------------ #
    plot_tools_panel(Xmax, Xaxis, df, CONTIG_ID, X, len)




    # ------------------------------ plot vlp coverage breadth ------------------------------ #
    vlp_cov_breadth_panel(Xmax, Xaxis, df, CONTIG_ID, W, B, study_col)




    # ------------------------------ plot vlp coverage depth ------------------------------ #
    vlp_cov_depth_panel(Xmax, Xaxis, df, CONTIG_ID, W, D, study_col, 4)

    legend(
        lty = 1,
        col = sapply(names(study_col), function (x) study_col[ x ]),
        legend = sub('VanEspen', 'Van Espen', sub('_', ' et al. ', names(study_col))),
        title = 'VLP reads:',
        'topright', inset = c(-0.01, -0.15), xpd = T, ncol = 1
    )




    # ------------------------------ plot own reads coverage depth ------------------------------ #
    own_reads_panel(Xmax, Xaxis, df, CONTIG_ID, Wq, Q, Q_range)

    legend(
        lty = 1,
        col = 'gray25',
        legend = 'own reads',
        'topright', inset = c(-0.01, -0.15), xpd = T
    )

    mtext(paste0(CONTIG_ID_short, ', kb'), side = 1, line = 3, at = max(Xaxis) / 2)

}

dev.off()
