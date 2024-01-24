dp <- function (file, yID, xleft = 0, ybottom = 0, col = 'black') {

    t <- tryCatch(
        read.table(file, sep = '\t', header = F, stringsAsFactors = F),
        error = function (e) { return(data.frame(NULL)) }
    )

    if (nrow(t) > 0) {

        x_rows <- which(t[, 1] != yID)
        y_rows <- which(t[, 1] == yID)

        d <- data.frame(
            x_from = t[x_rows, 4],
            x_to   = t[x_rows, 5],
            y_from = t[y_rows, 4],
            y_to   = t[y_rows, 5],
            stringsAsFactors = F
        )  

        apply(d, 1, function (v) lines(x = xleft + v[1:2], y = ybottom + v[3:4], lwd = 0.1, col = col))

    }

}
