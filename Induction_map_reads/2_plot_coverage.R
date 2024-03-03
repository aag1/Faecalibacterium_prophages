sessionInfo()



P <- read.table('../CPs_refine_coo/prophages_refined_coo.txt', sep = '\t', header = T, stringsAsFactors = F)

DF <- read.table('../Induction_raw_reads/Novogene_summary_table.txt', sep = '\t', header = T, stringsAsFactors = F)



Xmax <- 0

B <- list()

for (STRAIN_ID in unique(DF$Strain)) {

    f <- paste0('../Fprau_map_reads/BED_FILES/', STRAIN_ID, '_contigs.bed')

    b <- read.table(f, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)
    b <- b[b[, 2] > 50000, ]

    Xmax <- max(Xmax, sum(b[, 2]))

    B[[ STRAIN_ID ]] <- b

}



pdf('prophage_induction_mapped_reads.pdf', height = 7, width = 7.5)

options(scipen = 1000)

layout(matrix(1:6, nrow = 6))

par(mar = c(2.5, 4, 2, 1), oma = c(0, 1.25, 0, 0))

for (STRAIN_ID in names(B)) {

    ### contig lengths
    b <- B[[ STRAIN_ID ]]


    ### coverage depth
    samples <- DF$Sample.Name[(DF$Strain == STRAIN_ID) & (DF$Sequencing == 'Yes')]
    samples <- samples[c(1, length(samples):2)]

    L <- lapply(samples, function (x) {

        f <- paste0('MAPPED_READS/', x, '_mapped_to_', STRAIN_ID, '.depth.txt')

        read.table(f, sep = '\t', header = F, stringsAsFactors = F)

    })

    names(L) <- samples


    ### plot
    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(1, 10^5),
        log = 'y',
        xaxs = 'i',
        yaxs = 'i',
        bty = 'n',
        axes = F,
        ann = F
    )

    SI <- sub('_S[0-9]+$', '', STRAIN_ID)
    if (SI == 'FM7') { SI <- 'A2-165' }
    mtext(SI, side = 3, line = 0.25, at = 0, adj = 0, font = 2, cex = 0.75)

    axis(side = 2, las = 1, mgp = c(3, 0.5, 0), tck = -0.05, cex.axis = 0.75)

    mtext('contig', side = 1, line = 0.25, at = -25000, adj = 1, cex = 0.5)

    Xleft <- 0

    for (CONTIG_ID in rownames(b)) {  # for each contig ...

        ### prophages
        idx <- which((P$strain_id == STRAIN_ID) & (P$contig_id == CONTIG_ID))

        if (length(idx) > 0) {

            rect(
                xleft = Xleft + P$new_from[idx],
                xright = Xleft + P$new_to[idx],
                ybottom = 1,
                ytop = 10^5,
                col = 'grey85',
                border = NA,
                xpd = T
            )

            mtext(
                P$candidate_prophage[idx],
                side = 3,
                line = ifelse(P$candidate_prophage[idx] %in% c('CP5', 'CP6', 'CP7'), -1.25, 0.25),
                adj = ifelse(P$candidate_prophage[idx] == 'CP6', 0, 0.5),
                at = Xleft + sum(P$new_from[idx] + P$new_to[idx]) / 2,
                cex = 0.75
            )

        }


        ### coverage
        len <- b[CONTIG_ID, 2]

        for (SAMPLE_ID in samples) {    # for each sequencing sample ...

            COL <- 'blueviolet'
            idx <- which(DF$Sample.Name == SAMPLE_ID)
            if ((DF$DNA_type[idx] == 'bacterial DNA')) { COL <- 'cornflowerblue' }
            if ((DF$DNA_type[idx] == 'VLP DNA') & (DF$Treatment[idx] == 'none')) { COL <- 'firebrick3' }
            if (DF$Treatment[idx] == 'heat 40 C') { COL <- 'goldenrod1' }

            t <- L[[ SAMPLE_ID ]]

            i <- which(t[, 1] == CONTIG_ID)
            coo <- t[i, 2]
            cov <- t[i, 3]

            COVERAGE <- rep(0, len)
            COVERAGE[coo] <- cov

            w_center <- seq(from = 1501, to = len - 1500, by = 500)
            w_from <- w_center - 1500
            w_to <- w_center + 1500
            w_to[ length(w_to) ] <- len

            Q <- sapply(seq_along(w_center), function (i) { mean( COVERAGE[ w_from[i]:w_to[i] ] ) })

            lines(Xleft + w_center, Q + 1, col = COL)

        }

        rect(xleft = Xleft, xright = Xleft + len, ybottom = 1, ytop = 10^5)

        mtext(sub('^NODE_([0-9]+)_.+$', '\\1', CONTIG_ID), side = 1, line = 0.25, at = Xleft + len / 2, cex = 0.5)

        Xleft <- Xleft + len

    }


    ### scale
    if (STRAIN_ID == 'HTF-128_S3') {

        axis(
            side = 1,
            pos = 10^5,
            at = par()$usr[2] - c(150000, 50000),
            labels = c('0', '100 kb'),
            cex.axis = 0.75,
            mgp = c(3, 0.5, 0),
            xpd = T
        )

    }


    ### legend
    if (STRAIN_ID == 'FM8_S2') {

        legend(
            'right',
            legend = c('Bacterial DNA, no treatment', 'VLP DNA, no treatment', 'VLP DNA, 40 °C treatment', 'VLP DNA, MMC treatment'),
            col = c('cornflowerblue', 'firebrick3', 'goldenrod1', 'blueviolet'),
            lty = 1,
            lwd = 1.5,
            ncol = 1,
            cex = 0.75,
            xpd = T
        )

    }

}

par(new = T, oma = c(0, 1.25, 0, 0)); layout(1)

mtext('Mean depth + 1', side = 2, outer = T)

dev.off()
