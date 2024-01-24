# library TeachingDemos must be loaded!

plot_genome_map <- function (lab, len, tab, dom, dgr, dom_lab = T, dgr_col = 'black', cex_lab = 0.5, xLeft = 0, yBottom = 0) {

    # ------------------------------ ORFs ------------------------------ #
    D <- data.frame(NULL)


    for (i in 1:nrow(tab)) {

        from   <- tab[i, 2]
        to     <- tab[i, 3]
        strand <- tab[i, 4]


        if (strand == '+') {
            frame <- ((from - 1) %% 3) + 1
            Y <- 6 - frame
        }
        if (strand == '-') {
            frame <- ((len - to) %% 3) + 1
            Y <- 3 - frame
        }


        # grey ORF rectangle
        rect(
            xleft   = xLeft + from,
            xright  = xLeft + to,
            ybottom = yBottom + Y,
            ytop    = yBottom + Y + 1,
            col     = 'grey98',
            border  = NA
        )


        # protein domains
        idx <- which(dom$protein_id == tab[i, 1])

        if (length(idx) == 1) {

            s <- strsplit(dom$protein_env_coo[idx], ';')[[1]]
            l <- lapply(strsplit(s, '-'), as.numeric)
            c1 <- unlist(lapply(l, function (v) v[1]))
            c2 <- unlist(lapply(l, function (v) v[2]))

            if (strand == '+') {
                d_from <- (from - 1) + c1*3 - 2
                d_to <- (from - 1) + c2*3
            }
            if (strand == '-') {
                d_from <- (to + 1) - c2*3
                d_to <- (to + 1) - c1*3 + 2
            }

            rect(
                xleft   = xLeft + d_from,
                xright  = xLeft + d_to,
                ybottom = yBottom + Y,
                ytop    = yBottom + Y + 1,
                col     = dom$color[idx],
                border  = NA
            )

            d <- data.frame(
                lab = dom$profile_name[idx],
                coo = xLeft + (from + to) / 2,
                col = dom$color[idx],
                stringsAsFactors = F
            )
            D <- rbind(D, d)

        }


        # black ORF border
        rect(
            xleft   = xLeft + from,
            xright  = xLeft + to,
            ybottom = yBottom + Y,
            ytop    = yBottom + Y + 1,
            col     = NA,
            border  = 'black'
        )

    }


    # protein domain labels
    if (dom_lab) {

        D$coo <- spread.labs(D$coo, mindiff = 1000, maxiter = 10^6)

        text(
            D$lab,
            x = D$coo,
            y = yBottom + 7,
            col = D$col,
            cex = cex_lab,
            srt = 45,
            adj = 0,
            xpd = T
        )

    }



    # ------------------------------ DGR repeats ------------------------------ #
    if (!is.null(dgr)) {

        for (i in 1:nrow(dgr)) {

            TR_center <- sum(dgr[i, c('TR_coo1', 'TR_coo2')]) / 2
            VR_center <- sum(dgr[i, c('VR_coo1', 'VR_coo2')]) / 2

            lines(
                x = xLeft + rep(TR_center, 2),
                y = c(yBottom, yBottom - 2),
                col = dgr_col,
                xpd = T
            )

            lines(
                x = xLeft + c(TR_center, VR_center),
                y = rep(yBottom - 2, 2),
                col = dgr_col,
                xpd = T
            )

            arrows(
                x0 = xLeft + VR_center,
                y0 = yBottom - 2,
                y1 = yBottom,
                length = 0.08,
                col = dgr_col,
                xpd = T
            )

        }

    }



    # ------------------------------ genome ------------------------------ #
    rect(
        xleft   = xLeft,
        xright  = xLeft + len, 
        ybottom = yBottom,
        ytop    = yBottom + 6,
        col     = NA,
        border  = 'black'
    )


    text(
        c('+1', '+2', '+3', '-1', '-2', '-3'),
        x = xLeft - diff(par()$usr[1:2]) * 0.01,
        y = yBottom + 6:1 - 0.5,
        cex = cex_lab,
        adj = 1,
        xpd = T
    )


    text(
        lab,
        x = xLeft - diff(par()$usr[1:2]) * 0.05,
        y = yBottom + 3,
        adj = 1,
        xpd = T
    )

}
