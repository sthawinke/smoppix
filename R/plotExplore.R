#' Plot the hyperframe for exploratory purposes
#'
#' @param hypFrame The hyperframa
#' @param features A small number of features to be fitted. Defaults to the first 5
#' @param ppps The rownames or indices of the point patterns to be plotted. Defaults to the first 9.
#'
#' @return Plots a facet of point patterns to output
#' @importFrom spatstat.geom is.hyperframe
#' @importFrom grDevices palette
#' @importFrom graphics par legend
#' @export
#'
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c("x", "y"),
#'     imageVars = c("day", "root", "section")
#' )
#' plotExplore(hypYang)
plotExplore <- function(hypFrame, features = attr(hypFrame, "features"), ppps) {
    stopifnot(is.hyperframe(hypFrame), is.character(features))
    npp <- nrow(hypFrame)
    if (missing(ppps)) {
        ppps <- rownames(hypFrame)[seq_len(min(8, npp))]
    }
    Cols <- palette()
    Cols <- c(Cols, rep("grey", length(features) - length(Cols)))
    names(Cols) <- features
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mfrow = if ((LL <- length(ppps)) <= 3) c(2, 2) else c(3, 3), mar = c(3, 3, 3, 3))
    baa <- lapply(ppps, function(i) {
        id <- marks(hypFrame$ppp[[i]], drop = FALSE)$gene %in% features
        if (any(id)) {
            colVec <- Cols[marks(hypFrame$ppp[[i]], drop = FALSE)$gene[id]]
            cordMat <- coords(hypFrame$ppp[[i]][id, ])
            ordVec <- order(colVec != "grey") # Plot grey background first and coloured dots on top
            plot(cordMat[ordVec, ],
                main = ppps[i], pch = ".", asp = 1,
                col = colVec[ordVec]
            )
        }
    })
    plot(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    idCols <- which(Cols != "grey")
    legend("center",
        legend = names(Cols)[idCols], pch = 20, col = Cols[idCols],
        ncol = 2
    )
}
