#' Plot the hyperframe for exploratory purposes
#' @description A subset of 8 point patterns is plotted. This function is mainly
#' meant for exploratory purposes and confirming correct data read-in.
#'
#' @param hypFrame The hyperframe
#' @param features A small number of features to be fitted. Defaults to the first 5
#' @param ppps The rownames or indices of the point patterns to be plotted.
#' Defaults to the first 8.
#' @param maxPlot The maximum number of events plotted per point pattern
#' @param Cex Point amplification factor
#'
#' @return Plots a facet of point patterns to output
#' @importFrom spatstat.geom is.hyperframe coords
#' @importFrom grDevices palette
#' @importFrom graphics par legend
#' @export
#'
#' @examples
#' data(Yang)
#' hypYang <- buildHyperFrame(Yang,
#'     coordVars = c('x', 'y'),
#'     imageVars = c('day', 'root', 'section')
#' )
#' plotExplore(hypYang)
#' plotExplore(hypYang, features = c("SmRBRb", "SmTMO5b", "SmWER--SmAHK4f"))
plotExplore <- function(hypFrame, features = getFeatures(hypFrame)[seq_len(6)], ppps,
                        maxPlot = 1e+05, Cex = 1) {
    if(!is.hyperframe(hypFrame)){
        hypFrame <- hypFrame$hypFrame
    }
    stopifnot(is.hyperframe(hypFrame), is.character(features), is.numeric(maxPlot))
    features <- unique(unlist(lapply(features, sund)))
    npp <- nrow(hypFrame)
    if (missing(ppps)) {
        ppps <- seq_len(min(29, npp))
    }
    Cols <- palette()
    if(length(Cols) > length(features)){
        Cols <- Cols[seq_along(features)]
    }
    names(Cols) = features
    allFeats = getFeatures(hypFrame)
    if (length(Cols) < (lf <- length(allFeats))){
        Cols <- c(Cols, rep("grey", lf - length(Cols)))
        names(Cols)[Cols == "grey"] <- setdiff(allFeats, features)
    }
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mfrow = if ((LL <- length(ppps)) <= 3)
        c(2, 2) else if(LL <=8) c(3, 3) else if(LL <=15)
            c(4, 4) else if(LL <= 24) c(5, 5) else if(LL <= 29)
                c(5, 6) else if(LL <= 35) c(6,6) else c(6,7), mar = c(3, 3, 3, 3))
    baa <- lapply(ppps, function(i) {
            colVec <- Cols[marks(hypFrame$ppp[[i]], drop = FALSE)$gene]
            cordMat <- coords(hypFrame$ppp[[i]])
            ordVec <- order(colVec != "grey")
            # Plot grey background first and coloured dots on top
            plot(cordMat[ordVec, ], main = hypFrame$image[ppps[i]], pch = ".",
                 asp = 1, col = colVec[ordVec], cex = Cex)
    })
    plot(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    idCols <- which(Cols != "grey")
    legend("center", legend = names(Cols)[idCols], pch = 20, col = Cols[idCols], ncol = 2)
}
