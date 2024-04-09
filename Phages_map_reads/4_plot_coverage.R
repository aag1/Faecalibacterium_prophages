source('../Phages_genome_maps/function_plot_genome_map.R')
sessionInfo()




P <- read.table('../Induction_assembly/induced_phages_2.txt', sep = '\t', header = T, stringsAsFactors = F)

DF <- read.table('../Induction_raw_reads/Novogene_summary_table.txt', sep = '\t', header = T, stringsAsFactors = F)




pdf('phages_read_coverage.pdf', height = 8.75, width = 7.5)

layout(matrix(1:10, nrow = 10))

for (phage in c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie')) {

    # ------------------------------ load genome data ------------------------------ #
    load(paste0('../Phages_genome_maps/', phage, '_genome_map.RData'))



    # ------------------------------ coverage depth ------------------------------ #
    strain <- P$strain[ P$phage_name == phage ][1]
    samples <- DF$Sample.Name[(DF$Strain == strain) & (DF$DNA_type == 'VLP DNA') & (DF$Sequencing == 'Yes')]
    samples <- samples[c(length(samples):1)]

    L <- lapply(samples, function (x) {

        f <- paste0('MAPPED_READS/', x, '_mapped_to_', phage, '.depth.txt')

        read.table(f, sep = '\t', header = F, stringsAsFactors = F)

    })

    names(L) <- samples



    # ------------------------------ plot coverage depth ------------------------------ #
    par(mar = c(0, 6, 2, 1))

    if (phage == 'Tulp') { Yrange <- c(0, 1500) }
    if (phage == 'Roos') { Yrange <- c(10000, 60000) }
    if (phage == 'Pioen') { Yrange <- c(0, 12000) }
    if (phage == 'Aster') { Yrange <- c(0, 500) }
    if (phage == 'Lelie') { Yrange <- c(0, 1200) }

    w_center <- seq(from = 51, to = len - 50, by = 20)
    w_from <- w_center - 50
    w_to <- w_center + 50
    w_to[ length(w_to) ] <- len

    plot(
        NA,
        xlim = c(0, 55000),
        ylim = Yrange,
        xaxs = 'i', xaxt = 'n',
        yaxp = c(Yrange, ifelse(phage %in% c('Roos', 'Aster'), 5, 3)),
        las = 1, mgp = c(3.5, 0.8, 0.2), tck = -0.05, cex.axis = 0.8,
        ylab = 'Mean depth',
        bty = 'n',
    )

    for (SAMPLE_ID in samples) {

        COL1 <- 'blueviolet'
        idx <- which(DF$Sample.Name == SAMPLE_ID)
        if (DF$Treatment[idx] == 'none') { COL1 <- 'firebrick3' }
        if (DF$Treatment[idx] == 'heat 40 C') { COL1 <- 'goldenrod1' }

        t <- L[[ SAMPLE_ID ]]

        coo <- t[, 2]
        cov <- t[, 3]

        COVERAGE <- rep(0, len)
        COVERAGE[coo] <- cov

        Q <- sapply(seq_along(w_center), function (i) { mean( COVERAGE[ w_from[i]:w_to[i] ] ) })

        lines(w_center, Q + 1, col = COL1)

    }



    # ------------------------------ legend ------------------------------ #
    if (phage == 'Tulp') {

        legend(
            'topright',
            legend = c('VLP DNA, no treatment', 'VLP DNA, 40 °C treatment', 'VLP DNA, MMC treatment'),
            col = c('firebrick3', 'goldenrod1', 'blueviolet'),
            lty = 1,
            lwd = 1.5,
            ncol = 1,
            xpd = T
        )

    }



    # ------------------------------ plot genome map ------------------------------ #
    par(mar = c(3, 6, 1, 1))

    plot(NA, xlim = c(0, 55000), ylim = c(0, 6), xaxs = 'i', ann = F, axes = F)

    if (lab == 'Tulp') { lab <- 'Mushu' }
    plot_genome_map(lab, len, tab, dom, dgr = NULL, dom_lab = F, yBottom = 0)



    # ------------------------------ pac/cos site ------------------------------ #
    if (phage != 'Tulp') {

        points(x = 0, y = -1.5, pch = 17, cex = 1.25, xpd = T)

        text(ifelse(phage %in% c('Roos', 'Pioen'), 'pac', 'cos'), x = 0, y = -3.5, font = 3, xpd = T)

    }



    # ------------------------------ scale ------------------------------ #
    if (phage == 'Tulp') {

        axis(
            side = 1,
            pos = par()$usr[4],
            at = c(44000, 54000),
            labels = c('0', '10 kb')
        )

    }

}

dev.off()
