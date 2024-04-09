sessionInfo()




# ------------------------------ read data ------------------------------ #
tab <- data.frame(NULL)

for (p in c('Tulp', 'Roos', 'Pioen', 'Lelie', 'CP11', 'CP2')) {

    t <- read.table(paste0('../DGR_detection/', p, '_DGR_repeats.txt'), sep = '\t', header = T, stringsAsFactors = F)

    tab <- rbind(tab, t)

}

tab <- tab[, c('VR_orf', 'VR_aa1', 'VR_aa2', 'VR_aaL')]

tab$source <- 'phage'

tab$name <- sub('_([0-9]+)$', ' orf\\1', tab$VR_orf)
tab$name <- sub('^Tulp', 'Mushu', tab$name)
tab$name <- sub('^CP2', 'Lagaffe_CP', tab$name)


tab1 <- data.frame(
    VR_orf = c('DUF6273', 'Big_3_3', '6HHK_A', '1YU0_A'),
    VR_aa1 = c(NA, NA, NA, 337),
    VR_aa2 = c(NA, NA, NA, 381),
    VR_aaL = c(168, 157, 193, 381),
    source = c('pfam', 'pfam', 'pdb', 'pdb'),
    stringsAsFactors = F
)

tab1$name <- paste(c('', 'Ig-like', 'Tail fiber', 'MTD'), tab1$VR_orf)

tab <- rbind(tab, tab1)




# ------------------------------ palette ------------------------------ #
colfunc <- colorRampPalette(c('grey90', 'black'))
pal <- colfunc(101)




# ------------------------------ plot ------------------------------ #
pdf('DGR_targets_dotplots.pdf', height = 8.75, width = 7.5)

par(mar = c(7.5, 8, 0, 0), las = 2, tck = -0.01, cex.axis = 0.55, mgp = c(3, 0.5, 0))

plot(
    NA,
    xlim = c(0, sum(tab$VR_aaL[tab$source == 'phage']) + (sum(tab$source == 'phage') - 1) * 150),
    ylim = c(0, sum(tab$VR_aaL) + (nrow(tab) - 1) * 150 + 300),
    asp = 1, ann = F, axes = F
)


XL <- 0

for (i in rev(which(tab$source == 'phage'))) {

    XR <- XL + tab$VR_aaL[i]


    ### X axis
    text(tab$name[i], x = (XL + XR) / 2, y = -400, cex = 0.9, srt = 90, adj = 1, xpd = T)

    axis(side = 1, pos = 0, at = c(XL, XR), labels = c(0, tab$VR_aaL[i]))


    YB <- 0

    for (j in nrow(tab):1) {

        YT <- YB + tab$VR_aaL[j]


        ### VRs
        rect(xleft = XL + tab$VR_aa1[i], xright = XL + tab$VR_aa2[i], ybottom = YB, ytop = YT, col = 'papayawhip', border = NA)

        if (!is.na(tab$VR_aa1[j])) {

            rect(xleft = XL, xright = XR, ybottom = YB + tab$VR_aa1[j], ytop = YB + tab$VR_aa2[j], col = 'papayawhip', border = NA)

        }


        ### Y axis
        if (i == 1) {

            text(tab$name[j], x = -400, y = (YB + YT) / 2, cex = 0.9, adj = 1, xpd = T)

            axis(side = 2, pos = 0, at = c(YB, YT), labels = c(0, tab$VR_aaL[j]))

        }


        ### dotplot
        base <- paste0('HHalign_output/', tab$VR_orf[i], '___', tab$VR_orf[j])

        P <- read.table(paste0(base, '.probab'), header = F, stringsAsFactors = F)[, 1]
        P <- as.numeric(P)

        K <- read.table(paste0(base, '.atab1'), sep = '\t', header = F, stringsAsFactors = F)
        delim <- which(K[, 1] == 'i')
        if (length(P) != length(delim)) { stop('Different number of alignments in hhr and atab files!') }

        for (q in seq_along(delim)) {

            if (delim[q] == nrow(K)) { break }
            if ((q < length(delim)) & (delim[q] + 1 == delim[q + 1])) { next }

            r1 <- delim[q] + 1
            r2 <- ifelse(q == length(delim), nrow(K), delim[q + 1] - 1)

            k <- K[r1:r2, ]

            lines(
                x = XL + as.numeric(k[, 1]),
                y = YB + as.numeric(k[, 2]),
                col = pal[round(P[q]) + 1],
                lwd = 2
            )

        }


        ### frame
        rect(xleft = XL, xright = XR, ybottom = YB, ytop = YT)


        YB <- YT + ifelse((tab$source[j] != 'phage') & (tab$source[j - 1] == 'phage'), 450, 150)

    }


    XL <- XR + 150

}




# ------------------------------ legend ------------------------------ #
mtext(c('A', 'B'), side = 2, line = 6, at = c(7000, 1500), cex = 2.5)

for (i in 1:101) {

    rect(
        xleft  = par()$usr[1] - 1000 + (i - 1) * 7.5,
        xright = par()$usr[1] - 1000 + i * 7.5,
        ybottom = -800,
        ytop = -500,
        col = pal[i],
        border = NA,
        xpd = T
    )

}

text('Probability', x = par()$usr[1] - 1000 + 101 * 7.5 / 2, y = -400, cex = 0.9, xpd = T)

text(c(0, 100), x = par()$usr[1] - 1000 + c(0, 101 * 7.5), y = rep(-900, 2), cex = 0.55, xpd = T)

dev.off()
