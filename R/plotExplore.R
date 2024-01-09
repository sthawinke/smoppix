#' Plot the hyperframe for exploratory purposes
#' @description A subset of 8 point patterns is plotted. This function is mainly
#' meant for exploratory purposes and confirming correct data read-in.
#' @details The colour palette is taken from the output of palette(),
#' so set that one to change the colour scheme
#' @param hypFrame The hyperframe
#' @param features A small number of features to be fitted. Defaults to the first 5
#' @param ppps The rownames or indices of the point patterns to be plotted.
#' Defaults to the first 8.
#' @param maxPlot The maximum number of events plotted per point pattern
#' @param Cex Point amplification factor
#' @param plotWindows A boolean, should windows be plotted too?
#' @param Xlim,Ylim Vectors of length 2 with plotting limits
#'
#' @return Plots a facet of point patterns to output
#' @importFrom spatstat.geom is.hyperframe coords plot.owin
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
                        maxPlot = 1e+05, Cex = 1, plotWindows = !is.null(hypFrame$owins),
                        Xlim = NULL, Ylim = NULL) {
    if(!is.hyperframe(hypFrame)){
        hypFrame <- hypFrame$hypFrame
    }
    if(plotWindows && is.null(hypFrame$owins)){
        warning("No windows present in hyperframe object")
    }
    stopifnot(is.hyperframe(hypFrame), is.character(features), is.numeric(maxPlot))
    features <- unique(unlist(lapply(features, sund)))
    npp <- nrow(hypFrame)
    if (missing(ppps)) {
        ppps <- seq_len(min(29, npp))
    }
    Cols <- makeCols(features, hypFrame)
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mfrow = if ((LL <- length(ppps)) == 1) c(1,2) else if ((LL <- length(ppps)) <= 3)
        c(2, 2) else if(LL <=8) c(3, 3) else if(LL <=11) c(3,4) else if(LL <=15)
            c(4, 4) else if(LL <= 24) c(5, 5) else if(LL <= 29)
                c(5, 6) else if(LL <= 35) c(6,6) else c(6,7), mar = c(3, 3, 3, 3))
    baa <- lapply(ppps, function(i) {
            colVec <- Cols[marks(hypFrame$ppp[[i]], drop = FALSE)$gene]
            cordMat <- coords(hypFrame$ppp[[i]])
            ordVec <- order(colVec != "grey")
            # Plot grey background first and coloured dots on top
            plot(cordMat[ordVec, ], main = hypFrame$image[ppps[i]], pch = ".",
                 asp = 1, col = colVec[ordVec], cex = Cex, xlim = Xlim, ylim = Ylim)
            if(plotWindows)
                foo <- lapply(hypFrame$owins[[i]], plot.owin, add = TRUE)
    })
    plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    addLegend(Cols)
}
addLegend = function(Cols, Shift = c(0,0)){
    idCols <- which(Cols != "grey");li <- length(idCols)
    points(x = rep(0.1, li) + Shift[1], SeqY <- seq(0.95, 0.05, length.out = li)+ Shift[2], pch = 20,
           col = Cols[idCols])
    text(x = rep(0.4, li) + Shift[1], SeqY, names(Cols)[idCols], adj = 1)
}
makeCols = function(features, hypFrame){
    Cols <- setdiff(palette(), "black")
    if(length(Cols) > length(features)){
        Cols <- Cols[seq_along(features)]
    }
    names(Cols) = features
    allFeats = getFeatures(hypFrame)
    if (length(Cols) < (lf <- length(allFeats))){
        Cols <- c(Cols, rep("grey", lf - length(Cols)))
        names(Cols)[Cols == "grey"] <- setdiff(allFeats, features)
    }
    return(Cols)
}
