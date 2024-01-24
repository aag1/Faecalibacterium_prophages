legend_genome_map <- function (xLeft, yBottom, COL, cex_lab = 1) {

    ### frame
    rect(
        xleft = xLeft,
        xright = xLeft + 10000,
        ybottom = yBottom,
        ytop = yBottom + 18,
        xpd = T
    )



    ### virion
    rect(
        xleft = xLeft + 1000,
        xright = xLeft + 3000,
        ybottom = yBottom + 14.5,
        ytop = yBottom + 15.5,
        col = COL['virion'],
        border = NA,
        xpd = T
    )

    text('virion', x = xLeft + 4000, y = yBottom + 15, adj = 0, cex = cex_lab, xpd = T)



    ### DGR RT
    rect(
        xleft = xLeft + 1000,
        xright = xLeft + 3000,
        ybottom = yBottom + 10.5,
        ytop = yBottom + 11.5,
        col = COL['RT'],
        border = NA,
        xpd = T
    )

    text('RT', x = xLeft + 4000, y = yBottom + 11, adj = 0, cex = cex_lab, xpd = T)



    ### DGR repeats
    lines(
        x = rep(xLeft + 3000, 2),
        y = yBottom + c(2, 4),
        col = COL['RT'],
        xpd = T
    )

    lines(
        x = xLeft + c(3000, 7000),
        y = rep(yBottom + 2, 2),
        col = COL['RT'],
        xpd = T
    )

    arrows(
        x0 = xLeft + 7000,
        y0 = yBottom + 2,
        y1 = yBottom + 4,
        length = 0.08,
        col = COL['RT'],
        xpd = T
    )

    text(
        c('TR', 'VR'),
        x = xLeft + c(3000, 7000),
        y = yBottom + 6,
        cex = cex_lab,
        col = COL['RT'],
        xpd = T
    )

}
