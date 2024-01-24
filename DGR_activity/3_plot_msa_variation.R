library(bio3d)
source('../Phages_genome_maps/function_plot_genome_map.R')
source('../Phages_genome_maps/function_legend_genome_map.R')
sessionInfo()




pdf(paste0('DGR_msa_variation.pdf'), height = 8.75, width = 7.5)

layout(matrix(1:9, ncol = 1))




for (x in c('Roos', 'CP2', 'CP11')) {

    # ------------------------------ subset MSA ------------------------------ #
    f <- paste0(x, '_cognate_genomes_msa.fasta')
    MSA <- read.fasta(f)$ali
    MSA <- toupper(MSA)

    idx <- which(MSA[x, ] != '-')
    MSA <- MSA[, idx]



    # ------------------------------ sliding window ------------------------------ #
    w_center <- seq(from = 51, to = ncol(MSA) - 50, by = 20)
    w_from <- w_center - 50
    w_to <- w_center + 50
    w_to[ length(w_to) ] <- ncol(MSA)
    w_len <- w_to - w_from + 1



    # ------------------------------ MSA columns divergence ------------------------------ #
    diverg <- apply(MSA, 2, function (v) {

        q <- table(v)
        s <- names(q)[q == max(q)][1]
        sum(v != s) / length(v)

    })

    w_diverg <- sapply(seq_along(w_center), function (i) {

        mean( diverg[ w_from[i]:w_to[i] ] )

    })



    # ------------------------------ highly conserved MSA columns ------------------------------ #
    S <- c('A', 'T', 'G', 'C', '-')

    conserv <- list()
    for (s in S) {

        nrow_msa <- nrow(MSA)

        conserv[[s]] <- apply(MSA, 2, function (v) { sum(v == s) / nrow_msa >= 0.9 })

    }

    w_conserv_A <- sapply(seq_along(w_center), function (i) {

        w_coo <- w_from[i]:w_to[i]

        n_A <- sum(conserv[['A']][w_coo])

        n_all <- sum(sapply(S, function (s) sum(conserv[[s]][w_coo])))

        n_A / n_all * 100

    })



    # ------------------------------ load genome data ------------------------------ #
    load(paste0('../Phages_genome_maps/', x, '_genome_map.RData'))



    # ------------------------------ panel A ------------------------------ #
    par(mar = c(0.1, 6.5, 2, 1), mgp = c(3.2, 1.2, 0.2))

    yMax <- max(w_diverg)
    yMax <- ceiling(yMax / 0.2) * 0.2

    plot(
        NA,
        xlim = c(0, 55000),
        ylim = c(0, yMax),
        xaxs = 'i',
        yaxs = 'i',
        xaxt = 'n',
        ylab = 'MSA\nvariation',
        bty = 'n',
        las = 1
    )

    for (i in 1:nrow(dgr)) {

        lines(
            x = rep((dgr$VR_coo1[i] + dgr$VR_coo2[i]) / 2, 2),
            y = c(-1, par()$usr[4]),
            col = 'papayawhip',
            lend = 'butt', lwd = 2,
            xpd = T
        )

    }

    lines(x = w_center, y = w_diverg, col = 'grey25')



    # ------------------------------ panel B ------------------------------ #
    yMax <- max(w_conserv_A, na.rm = T)
    yMax <- ceiling(yMax / 20) * 20

    plot(
        NA,
        xlim = c(0, 55000),
        ylim = c(0, yMax),
        xaxs = 'i',
        yaxs = 'i',
        xaxt = 'n',
        ylab = 'Conserved\nA, %',
        bty = 'n',
        las = 1
    )

    for (i in 1:nrow(dgr)) {

        lines(
            x = rep((dgr$VR_coo1[i] + dgr$VR_coo2[i]) / 2, 2),
            y = c(-100, 200),
            col = 'papayawhip',
            lend = 'butt', lwd = 2,
            xpd = T
        )

    }

    lines(x = w_center, y = w_conserv_A, col = 'grey25')



    # ------------------------------ panel C ------------------------------ #
    par(mar = c(3, 6.5, 2, 1))

    plot(
        NA,
        xlim = c(0, 55000),
        ylim = c(0, 6),
        xaxs = 'i',
        yaxs = 'i',
        axes = F,
        ann = F
    )

    for (i in 1:nrow(dgr)) {

        lines(
            x = rep((dgr$VR_coo1[i] + dgr$VR_coo2[i]) / 2, 2),
            y = c(6, 100),
            col = 'papayawhip',
            lend = 'butt', lwd = 2,
            xpd = T
        )

    }

    plot_genome_map(lab, len, tab, dom, dgr, dom_lab = F, dgr_col = COL['RT'], yBottom = 0)

}




axis(
    side = 1,
    pos = par()$usr[4],
    at = c(44000, 54000),
    labels = c('0', '10 kb'),
    mgp = c(3, 1, 0)
)

layout(1)
par(new = T, mar = c(5/3, 13/3, 40.2, 2/3))
plot(NA, xlim = c(0, 55000), ylim = c(0, 6), xaxs = 'i', yaxs = 'i', axes = F, ann = F)
legend_genome_map(44000, 15, COL, 2/3)

dev.off()
