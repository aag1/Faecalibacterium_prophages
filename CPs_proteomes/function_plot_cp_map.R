# library TeachingDemos must be loaded!

plot_cp_map <- function (lab, len, tab, dom, xLeft = 0, yBottom = 0) {

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
    D$coo <- spread.labs(D$coo, mindiff = 1000, maxiter = 10^6)

    text(
        D$lab,
        x = D$coo,
        y = yBottom + 6.5,
        col = D$col,
        cex = 0.25,
        srt = 45,
        adj = 0,
        xpd = T
    )
    



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
        cex = 0.25,
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
