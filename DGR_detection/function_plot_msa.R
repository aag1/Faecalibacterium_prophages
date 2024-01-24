plot_msa <- function (MSA, COL, xLeft = 0, yBottom = 0) {

    for (i in 1:nrow(MSA)) {

        sapply(1:ncol(MSA), function (j) {

            rect(
                xleft = xLeft + (j - 1),
                xright = xLeft + j,
                ybottom = yBottom + nrow(MSA) - i,
                ytop = yBottom + nrow(MSA) - i + 1,
                col = COL[ MSA[i, j] ],
                border = NA
            )

        })


        rect(
            xleft = xLeft,
            xright = xLeft + ncol(MSA),
            ybottom = yBottom + nrow(MSA) - i,
            ytop = yBottom + nrow(MSA) - i + 1
        )


        text(
            rownames(MSA)[i],
            x = xLeft - 2,
            y = yBottom + nrow(MSA) - i + 0.5,
            adj = 1,
            xpd = T
        )

    }


    axis(
        side = 3,
        pos = yBottom + nrow(MSA),
        at = seq(from = 0, to = floor(ncol(MSA) / 50) * 50, by = 50)
    )

}
