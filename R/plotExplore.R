#' Plot the hyperframe for exploratory purposes
#' @description A subset of point patterns is plotted that fit in the window, highlighting some features.
#' This function is meant for exploratory purposes as well as for visual confirmation of findings.
#' @note The colour palette is taken from the output of palette(),
#' so set that one to change the colour scheme.
#'
#' @param hypFrame The hyperframe
#' @param features A small number of features to be fitted. Defaults to the first 5
#' @param ppps The rownames or indices of the point patterns to be plotted.
#' Defaults to maximum 99.
#' @param maxPlot The maximum number of events plotted per point pattern
#' @param Cex Point amplification factor
#' @param plotWindows A boolean, should windows be plotted too?
#' @param plotPoints A boolean, should the molecules be plotted as points?
#' @param Cex.main Expansion factor for the title
#' @param Mar the margins
#' @param titleVar Image variable to be added to the title
#' @param Xlim,Ylim plotting limits
#' @param piColourCell PI by which to colour the cell
#' @param piEsts Set of PI estimates, returned by estPis
#' @param palCols Two extremes of the colour palette for colouring the cells
#' @param border Passed on to plot.owin, and further to graphics::polygon
#' @note palCols sets the pseudo-continuous scale to colour cells.
#'
#' @return Plots a facet of point patterns to output
#' @importFrom spatstat.geom is.hyperframe coords plot.owin
#' @importFrom grDevices palette
#' @importFrom graphics par
#' @importFrom stats aggregate
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @examples
#' example(buildHyperFrame, 'smoppix')
#' plotExplore(hypYang)
#' plotExplore(hypYang, titleVar = 'day')
#' plotExplore(hypYang, features = c('SmRBRb', 'SmTMO5b', 'SmWER--SmAHK4f'))
plotExplore <- function(hypFrame, features = getFeatures(hypFrame)[seq_len(6)], ppps,
    maxPlot = 1e+05, Cex = 1, plotWindows = !is.null(hypFrame$owins), plotPoints = TRUE, piEsts = NULL,
    Xlim = NULL, Ylim = NULL, Cex.main = 0.8, Mar = c(0.35, 0.1, 0.75, 0.1), titleVar = NULL,
    piColourCell = NULL, palCols = c("blue", "yellow"), border = NULL) {
    if (!is.hyperframe(hypFrame)) {
        hypFrame <- hypFrame$hypFrame
    }
    if (plotWindows && is.null(hypFrame$owins)) {
        warning("No windows present in hyperframe object")
    }
    stopifnot(is.hyperframe(hypFrame), is.character(features), is.numeric(maxPlot),
        is.null(titleVar) || titleVar %in% getPPPvars(hypFrame))
    if (colourCells <- !is.null(piColourCell) && plotWindows) {
        featsplit <- sund(features)
        if (is.null(piEsts)) {
            stop("Supply list of PI estimates to lmms argument to colour cells")
        }
        piDf <- buildDataFrame(piEsts, gene = featsplit, pi = piColourCell)
        if (piColourCell %in% c("edge", "centroid")) {
            piDf <- aggregate(pi ~ cell +image, piDf, FUN = mean, na.rm = TRUE)
        }
        foo <- checkPi(piEsts, piColourCell)
        piColourCell <- match.arg(piColourCell, choices = c("edge", "centroid", "nnCell",
            "nnPairCell"))
        if (length(features) != 1) {
            stop("Supply single gene for colouring cells by univariate PI!")
        }
        pal <- colorRampPalette(palCols)
        Order <- findInterval(piDf$pi, sort(piDf$pi))
        Palette <- pal(sum(!is.na(piDf$pi)))[Order]
    }
    features <- unique(unlist(lapply(features, sund)))
    npp <- nrow(hypFrame)
    if (missing(ppps)) {
        ppps <- seq_len(min(99, npp))
    } else if (is.character(ppps)) {
        ppps <- match(ppps, rownames(hypFrame))
    }
    Cols <- makeCols(features, hypFrame)
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    LL <- length(ppps)
    Nrow <- ceiling(sqrt(LL))
    Ncol <- (LL%/%Nrow) + 1
    par(mfrow = c(Nrow, Ncol), mar = Mar)
    baa <- lapply(ppps, function(i) {
        PPPsub <- subSampleP(if (colourCells) {
            hypFrame$ppp[[i]][marks(hypFrame$ppp[[i]])$gene %in% features, ]
        } else {
            hypFrame$ppp[[i]]
        }, maxPlot)
        colVec <- Cols[marks(PPPsub, drop = FALSE)$gene]
        cordMat <- coords(PPPsub)
        ordVec <- order(colVec != "grey")
        plot(cordMat[ordVec, ], main = paste(hypFrame$image[i], if (!is.null(titleVar)) {
            hypFrame[[i, titleVar]]
        }), type = "n", cex.main = Cex.main, xaxt = "n", yaxt = "n", xaxs = "i",
            yaxs = "i", asp = 1, xlim = Xlim, ylim = Ylim)
        if (plotWindows) {
            if(colourCells){
                piDf2 = piDf[id <- (rownames(hypFrame)[i] == piDf$image),]
                Palette2 = Palette[id]
            }
            foo <- lapply(names(hypFrame$owins[[i]]), function(cell) {
                plot.owin(hypFrame[i, "owins", drop = TRUE][[cell]], add = TRUE,
                  col = if (colourCells && length(Col <- Palette2[cell == piDf2$cell])
                    && !is.na(Col)) {
                    Col
                  }, border = border)
            })
        }
        if(plotPoints){
            points(cordMat[ordVec, ], pch = ".", col = colVec[ordVec], cex = Cex)
        }

    })
    plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    colPlot <- if (colourCells) {
        tmp <- pal(6)
        names(tmp) <- signif(seq(min(piDf$pi, na.rm = TRUE), max(piDf$pi, na.rm = TRUE),
            length.out = 6), 2)
        tmp
    } else {
        Cols
    }
    addLegend(colPlot, Main = if (colourCells) {
        paste("pi", piColourCell)
    } else {
        "Feature"
    })
}
addLegend <- function(Cols, Shift = c(0, 0), Cex = 0.85, Pch = 20, Main = "") {
    idCols <- which(Cols != "grey")
    li <- length(idCols)
    maxY <- 0.85 - 0.1 * (mainId <- Main != "")
    points(x = rep(0.15, li) + Shift[1], SeqY <- seq(maxY, 0.1, length.out = li) +
        Shift[2], pch = Pch, col = Cols[idCols])
    text(x = rep(0.25, li) + Shift[1], SeqY, names(Cols)[idCols], adj = 0, cex = Cex)
    if (mainId) {
        text(x = 0.4 + Shift[1], y = 0.95 + Shift[2], Main, cex = 1.2)
    }
}
makeCols <- function(features, hypFrame) {
    Cols <- setdiff(c("red", "blue", rev(palette("Set1"))), "black")
    if (length(Cols) > length(features)) {
        Cols <- Cols[seq_along(features)]
    }
    names(Cols) <- features
    allFeats <- getFeatures(hypFrame)
    if (length(Cols) < (lf <- length(allFeats))) {
        Cols <- c(Cols, rep("grey", lf - length(Cols)))
        names(Cols)[Cols == "grey"] <- setdiff(allFeats, features)
    }
    return(Cols)
}
