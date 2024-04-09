source('../Phages_genome_maps/function_plot_genome_map.R')
sessionInfo()




P <- read.table('../Induction_assembly/induced_phages_2.txt', sep = '\t', header = T, stringsAsFactors = F)

DF <- read.table('../Induction_raw_reads/Novogene_summary_table.txt', sep = '\t', header = T, stringsAsFactors = F)




plot_var <- function (x, SAMPLE_ID, strain, DF, w_center, w_from, w_to, axis_side) {

    COL1 <- 'gray25'
    if (SAMPLE_ID != strain) {

        COL1 <- 'blueviolet'
        idx <- which(DF$Sample.Name == SAMPLE_ID)
        if (DF$DNA_type[idx] == 'bacterial DNA') { COL1 <- 'cornflowerblue' }
        if ((DF$DNA_type[idx] == 'VLP DNA') & (DF$Treatment[idx] == 'none')) { COL1 <- 'firebrick3' }
        if (DF$Treatment[idx] == 'heat 40 C') { COL1 <- 'goldenrod1' }

    }



    f <- paste0('../Phages_map_reads/NT_VARIATION/', x, '_', SAMPLE_ID, '_nt_variation.txt')
    if (file.exists(f)) {

        t <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)



        t$nt_var <- sapply(1:nrow(t), function (i) {

            q <- unlist(t[i, c('A', 'T', 'G', 'C')])
            i <- which(q == max(q))[1]
            o <- q[-i]
            sum(o) / sum(q)

        })



        Q <- sapply(seq_along(w_center), function (i) {

            mean( t$nt_var[ t$pos %in% w_from[i]:w_to[i] ] )

        })



        yMax <- ifelse(max(Q) > 0.05, 0.1, 0.05)

        plot(
            NA,
            xlim = c(0, 55000),
            ylim = c(0, yMax),
            xaxs = 'i',
            yaxs = 'i',
            axes = F,
            ann = F,
            bty = 'n'
        )

        axis(
            side = axis_side,
            pos = ifelse(axis_side == 2, -500, len + 500),
            at = c(0, yMax),
            labels = c('0', yMax),
            las = 1,
            mgp = c(3, 0.6, 0), tck = -0.2, cex.axis = 0.8
        )

        for (i in 1:nrow(dgr)) {

            lines(
                x = rep((dgr$VR_coo1[i] + dgr$VR_coo2[i]) / 2, 2),
                y = c(-1, 2),
                col = 'papayawhip',
                lend = 'butt', lwd = 2,
                xpd = T
            )

        }



        lines(w_center, Q, col = COL1, xpd = T)

    }

}




pdf(paste0('DGR_reads_variation.pdf'), height = 8.75, width = 7.5)

layout(matrix(1:24, ncol = 1), height = c(rep(0.5, 4), 1.5, 0.75, rep(0.5, 8), 1.5, 0.75, rep(0.5, 3), 1.5, 0.75, 0.75, 0.5, 1.5))

par(oma = c(0, 2.5, 2, 0))



for (x in c('Tulp', 'Roos', 'Pioen', 'Lelie')) {

    # ------------------------------ load genome data ------------------------------ #
    load(paste0('../Phages_genome_maps/', x, '_genome_map.RData'))



    # ------------------------------ list samples ------------------------------ #
    strain <- P$strain[ P$phage_name == x ][1]
    if (x == 'CP2') { strain <- 'FM7_S1' }
    if (x == 'CP11') { strain <- 'HTF-128_S3' }
    samples <- DF$Sample.Name[(DF$Strain == strain) & (DF$Sequencing == 'Yes')]



    # ------------------------------ plot variation ------------------------------ #
    par(mar = c(0.5, 5, 0.1, 3))


    w_center <- seq(from = 51, to = len - 50, by = 20)
    w_from <- w_center - 50
    w_to <- w_center + 50
    w_to[ length(w_to) ] <- len


    axis_side <- 2

    for (SAMPLE_ID in c(strain, samples)) {

        plot_var(x, SAMPLE_ID, strain, DF, w_center, w_from, w_to, axis_side)

        axis_side <- ifelse(axis_side == 2, 4, 2)

    }



    # ------------------------------ plot genome map ------------------------------ #
    par(mar = c(3, 5, 0.1, 3))

    plot(NA, xlim = c(0, 55000), ylim = c(0, 6), xaxs = 'i', ann = F, axes = F)

    if (lab == 'Tulp') { lab <- 'Mushu' }
    plot_genome_map(lab, len, tab, dom, dgr, dgr_col = COL['RT'], dom_lab = F, yBottom = 0)



    # ------------------------------ scale ------------------------------ #
    if (x == 'Tulp') {

        axis(
            side = 1,
            pos = par()$usr[4],
            at = c(44000, 54000),
            labels = c('0', '10 kb')
        )

    }

}



mtext('Reads variation', side = 2, line = 0.5, outer = T)

layout(1)
par(new = T, mar = c(0, 0, 0, 1))
plot(NA, xlim = 0:1, ylim = 0:1, axes = F, ann = F)
legend(
    'topright',
    legend = c('Bacterial genome sequencing', 'Bacterial DNA, no treatment', 'VLP DNA, no treatment', 'VLP DNA, 40 °C treatment', 'VLP DNA, MMC treatment'),
    col = c('gray25', 'cornflowerblue', 'firebrick3', 'goldenrod1', 'blueviolet'),
    lty = 1, lwd = 1.5,
    ncol = 1,
    cex = 0.6
)

dev.off()




pdf(paste0('DGR_reads_variation_CPs.pdf'), height = 2.6, width = 7.5)

layout(matrix(1:5, ncol = 1), height = c(rep(0.6, 2), 1.5, 0.6, 1.5))

par(oma = c(0, 2.5, 2, 0))



for (x in c('CP2', 'CP11')) {

    # ------------------------------ load genome data ------------------------------ #
    load(paste0('../Phages_genome_maps/', x, '_genome_map.RData'))



    # ------------------------------ list samples ------------------------------ #
    if (x == 'CP2') { strain <- 'FM7_S1' }
    if (x == 'CP11') { strain <- 'HTF-128_S3' }
    samples <- DF$Sample.Name[(DF$Strain == strain) & (DF$Sequencing == 'Yes')]



    # ------------------------------ plot variation ------------------------------ #
    par(mar = c(0.5, 7.5, 0.1, 3))


    w_center <- seq(from = 51, to = len - 50, by = 20)
    w_from <- w_center - 50
    w_to <- w_center + 50
    w_to[ length(w_to) ] <- len


    axis_side <- 2

    for (SAMPLE_ID in c(strain, samples)) {

        plot_var(x, SAMPLE_ID, strain, DF, w_center, w_from, w_to, axis_side)

        axis_side <- ifelse(axis_side == 2, 4, 2)

    }



    # ------------------------------ plot genome map ------------------------------ #
    par(mar = c(3, 7.5, 0.1, 3))

    plot(NA, xlim = c(0, 55000), ylim = c(0, 6), xaxs = 'i', ann = F, axes = F)

    if (lab == 'CP2') { lab <- 'Lagaffe_CP' }
    plot_genome_map(lab, len, tab, dom, dgr, dgr_col = COL['RT'], dom_lab = F, yBottom = 0)



    # ------------------------------ scale ------------------------------ #
    if (x == 'CP11') {

        axis(
            side = 1,
            pos = par()$usr[3],
            at = c(44000, 54000),
            labels = c('0', '10 kb')
        )

    }

}



mtext('Reads variation', side = 2, line = 0.5, outer = T, adj = 0.6)

layout(1)
par(new = T, mar = c(3, 0, 0, 1))
plot(NA, xlim = 0:1, ylim = 0:1, axes = F, ann = F)
legend(
    'bottomright',
    legend = c('Bacterial genome sequencing', 'Bacterial DNA, no treatment'),
    col = c('gray25', 'cornflowerblue'),
    lty = 1, lwd = 1.5,
    ncol = 1,
    cex = 0.6
)

dev.off()
