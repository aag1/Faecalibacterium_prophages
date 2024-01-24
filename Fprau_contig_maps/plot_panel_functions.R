plot_ORFs <- function (tab, len) {

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


        rect(
            xleft = from,
            xright = to,
            ybottom = Y,
            ytop = Y + 1,
            col = ifelse(strand == '+', 'lightcoral', 'royalblue1'),
            border = NA
        )

    }


    rect(
        xleft = 0,
        xright = len, 
        ybottom = 0,
        ytop = 6,
        col = NA,
        border = 'black'
    )


    text(
        x = -0.01 * diff(par()$usr[1:2]),
        y = 6:1 - 0.5,
        labels = c('+1', '+2', '+3', '-1', '-2', '-3'),
        cex = 0.75,
        adj = 1,
        xpd = TRUE
    )

}




# ---------------------------------------------------------------------------------------------------- #
plot_ORFs_panel <- function (Xmax, df, CONTIG_ID, len, cex_lab = 1) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(0, 6),
        xaxs = 'i', yaxs = 'i', 
        axes = FALSE, ann = FALSE, bty = 'n'
    )


    if (nrow(df) > 0) {

        for (i in 1:nrow(df)) {

            mtext(
                df$candidate_prophage[i],
                side = 3,
                line = 0.5,
                at = df$from[i] + (df$to[i] - df$from[i] + 1) / 2,
                col = 'black',
                cex = cex_lab,
                xpd = T
            )

            rect(
                xleft = df$from[i],
                xright = df$to[i],
                ybottom = -6,
                ytop = 6,
                col = 'grey85',
                border = NA,
                xpd = T
            )

        }

    }


    plot_ORFs(P[P$V1 == CONTIG_ID, ], len)

}




# ---------------------------------------------------------------------------------------------------- #
plot_nt_content_panel <- function (Xmax, Xaxis, df, CONTIG_ID, Wn, N, nt_col, x_axis_labs = T) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(10, 40),
        xaxs = 'i', yaxs = 'i',
        xlab = '', ylab = 'Nucleotide, %',
        axes = F,
        bty = 'n'
    )


    if (nrow(df) > 0) {

        for (i in 1:nrow(df)) {

            rect(
                xleft = df$from[i],
                xright = df$to[i],
                ybottom = -20,
                ytop = 60,
                col = 'grey85',
                border = NA,
                xpd = T
            )

        }

    }


    for (NT_LETTER in names(nt_col)) {

        lines(x = Wn[[ CONTIG_ID ]], y = N[[ CONTIG_ID ]][[ NT_LETTER ]], col = nt_col[ NT_LETTER ])

    }


    axis(side = 1, at = Xaxis, labels = if (x_axis_labs) { Xaxis / 1000 } else { F })

    axis(side = 2)

}




# ---------------------------------------------------------------------------------------------------- #
plot_tools_panel <- function (Xmax, Xaxis, df, CONTIG_ID, X, len, cex_lab = 1, lwd_signal = 8, x_axis_labs = T) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(0.5, 3.5),
        xaxs = 'i', yaxs = 'i',
        ann = F,
        axes = F,
        bty = 'n'
    )


    if (nrow(df) > 0) {

        for (i in 1:nrow(df)) {

            rect(
                xleft = df$from[i],
                xright = df$to[i],
                ybottom = par()$usr[3] - diff(par()$usr[3:4]),
                ytop = par()$usr[4] + diff(par()$usr[3:4]),
                col = 'grey85',
                border = NA,
                xpd = T
            )

        }

    }


    text(c('Cenote-Taker2', 'VirSorter2', 'geNomad'), x = par()$usr[1] - diff(par()$usr[1:2]) * 0.01, y = 1:3, adj = 1, cex = cex_lab, xpd = T)

    for (i in 1:3) {

        lines(x = c(1, len), y = rep(i, 2), col = 'grey60', lty = 2)
    
        t <- X[[ i ]]

        t <- t[t$contig_id == CONTIG_ID, ]

        if (nrow(t) > 0) {
    
            for (j in 1:nrow(t)) {

                lines(
                    x = c(t$from[j], t$to[j]),
                    y = rep(i, 2),
                    col = 'gray25',
                    lend = 'butt',
                    lwd = lwd_signal
                )

            }

        }

    }


    axis(side = 1, at = Xaxis, labels = if (x_axis_labs) { Xaxis / 1000 } else { F })

}




# ---------------------------------------------------------------------------------------------------- #
vlp_cov_breadth_panel <- function (Xmax, Xaxis, df, CONTIG_ID, W, B, study_col, x_axis_labs = T) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(0, 100),
        xaxs = 'i', yaxs = 'i',
        xlab = '', ylab = 'Breadth, %',
        axes = F,
        bty = 'n'
    )


    if (nrow(df) > 0) {

        for (i in 1:nrow(df)) {

            rect(
                xleft = df$from[i],
                xright = df$to[i],
                ybottom = par()$usr[3] - diff(par()$usr[3:4]),
                ytop = par()$usr[4] + diff(par()$usr[3:4]),
                col = 'grey85',
                border = NA,
                xpd = T
            )

        }

    }


    for (STUDY_NAME in names(study_col)) {

        invisible(lapply(B[[ CONTIG_ID ]][[ STUDY_NAME ]], function (v) {

            lines(x = W[[ CONTIG_ID ]], y = v, col = study_col[ STUDY_NAME ], xpd = T)

        }))

    }


    axis(side = 1, at = Xaxis, labels = if (x_axis_labs) { Xaxis / 1000 } else { F })

    axis(side = 2)

}




# ---------------------------------------------------------------------------------------------------- #
vlp_cov_depth_panel <- function (Xmax, Xaxis, df, CONTIG_ID, W, D, study_col, max_power, x_axis_labs = T) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = c(1, 10^max_power),
        xaxs = 'i', yaxs = 'i',
        log = 'y',
        xlab = '', ylab = 'Mean depth + 1',
        axes = F,
        bty = 'n'
    )


    if (nrow(df) > 0) {

        for (i in 1:nrow(df)) {

            rect(
                xleft = df$from[i],
                xright = df$to[i],
                ybottom = 10^(-10),
                ytop = 10^10,
                col = 'grey85',
                border = NA,
                xpd = T
            )

        }

    }


    for (STUDY_NAME in names(study_col)) {

        invisible(lapply(D[[ CONTIG_ID ]][[ STUDY_NAME ]], function (v) {

            lines(x = W[[ CONTIG_ID ]], y = v + 1, col = study_col[ STUDY_NAME ], xpd = T)

        }))

    }


    axis(side = 1, at = Xaxis, labels = if (x_axis_labs) { Xaxis / 1000 } else { F })

    axis(side = 2, at = 10^c(0:max_power))

}




# ---------------------------------------------------------------------------------------------------- #
own_reads_panel <- function (Xmax, Xaxis, df, CONTIG_ID, Wq, Q, Q_range, x_axis_labs = T) {

    plot(
        NA,
        xlim = c(0, Xmax),
        ylim = Q_range,
        xaxs = 'i', yaxs = 'i',
        xlab = '', ylab = 'Mean depth',
        axes = F,
        bty = 'n'
    )


    if (nrow(df) > 0) {

        for (i in 1:nrow(df)) {

            rect(
                xleft = df$from[i],
                xright = df$to[i],
                ybottom = par()$usr[3],
                ytop = 2 * Q_range[2],
                col = 'grey85',
                border = NA,
                xpd = T
            )

        }

    }


    lines(x = Wq[[ CONTIG_ID ]], y = Q[[ CONTIG_ID ]], col = 'gray25', xpd = T)


    axis(side = 1, at = Xaxis, labels = if (x_axis_labs) { Xaxis / 1000 } else { F })

    axis(side = 2)

}
