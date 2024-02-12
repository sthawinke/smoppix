#' Plot the n cells with highest abundance of a feature, breaking their spatial orientation.
#' @description After testing for within-cell patterns, it may be useful to
#' look at the cells with the most events for certain genes. These are plotted
#' here, but the spatial location of the cells is lost! The choice and ranking of cells
#' is one of decreasing gene (pair) expression.
#'
#' @param obj A hyperframe, or an object containing one
#' @param features The features to be plotted, a character vector
#' @param nCells An integer, the number of cells to be plotted
#' @param Cex The point expansion factor
#' @param borderColVar The variable to colour borders of the cell
#' @param borderCols Colour palette for the borders
#' @param warnPosition A boolean, should a warning be issued printed on the
#'  image that cells are not in their original location?
#' @param ... Additional arguments, currently ignored
#' @inheritParams plotExplore
#'
#' @return Plots cells with highest expression to the plotting window, returns invisible
#' @export
#' @importFrom spatstat.geom plot.owin subset.ppp
#' @importFrom graphics points text
#'
#' @examples
#' example(addCell, "spatrans")
#' plotCells(hypFrame2, "gene1")
#' plotCells(hypFrame2, "gene1", borderColVar = "condition")
plotCells <- function(obj, features = getFeatures(obj)[seq_len(3)],
                      nCells = 100, Cex = 1.5, borderColVar = NULL, borderCols = rev(palette()),
                      Mar = c(0.35, 0.1, 0.75, 0.1), warnPosition = TRUE, ...) {
    if (!is.hyperframe(obj)) {
        obj <- obj$hypFrame
    }
    stopifnot(is.hyperframe(obj), !is.null(obj$owins), length(nCells) == 1,
              all(features %in% getFeatures(obj)))
    old.par <- par(no.readonly = TRUE);on.exit(par(old.par))
    par(mar = Mar)
    features <- unique(unlist(lapply(features, sund)))
    Cols <- makeCols(features, obj)
    tablesCell <- lapply(obj$ppp, function(p) {
        table(marks(p[marks(p, drop = FALSE)$gene %in% features, ], drop = FALSE)$cell)
    })
    names(tablesCell) <- rownames(obj)
    nCells <- min(nCells - 1, length(ul <- unlist(tablesCell)))
    nthcell <- sort(ul, decreasing = TRUE)[min(nCells + 1, length(ul))]
    tablesCell <- lapply(tablesCell, function(x) {
        x[x > nthcell & names(x) != "NA"]
    })
    if (colourBorder <- !is.null(borderColVar)) {
        if (borderColVar %in% getEventVars(obj)) {
            unVals <- unique(unlist(lapply(obj$ppp, function(x)
                unique(marks(x, drop = FALSE)[[borderColVar]]))))
            names(borderCols)[seq_along(unVals)] <- unVals
            borderCols <- lapply(seq_along(tablesCell), function(x) {
                borderCols[marks(obj$ppp[[x]])[[borderColVar]][match(marks(obj$ppp[[x]])$cell, names(tablesCell[[x]]))]]
            })
        } else if (borderColVar %in% getPPPvars(obj)) {
            unVals <- unique(obj[[borderColVar]])
            names(borderCols)[seq_along(unVals)] <- unVals
            borderCols <- lapply(seq_along(tablesCell), function(x) {
                rep(borderCols[[obj[[x, borderColVar]]]], length(tablesCell[[x]]))
            })
        } else {
            stop("Border colour variable not found!")
        }
    } else {
        borderCols <- lapply(seq_along(tablesCell), function(x) {
            rep("black", length(tablesCell[[x]]))
        })
    }
    counter <- 0L
    Ceils = splitWindow(nCells + 1 + colourBorder)
    plot(type = "n", x = c(0, Ceils[1]), y = c(0, Ceils[2]), xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = FALSE)
    for (i in seq_along(tablesCell)) {
        nam <- names(tablesCell)[i]
        ppp <- subset.ppp(obj[[nam, "ppp"]], gene %in% features)
        for (j in seq_along(tablesCell[[nam]])) {
            namIn <- names(tablesCell[[nam]])[j]
            shifted <- toUnitSquare(obj[[nam, "owins"]][[namIn]], ppp = subset.ppp(ppp, cell == namIn), Shift = Shift <- c(
                counter %% Ceils[1],
                counter %/% Ceils[1]
            ))
            plot.owin(shifted$owin, add = TRUE, border = borderCols[[i]][[j]])
            # text(x = Shift[1], y = Shift[2], labels = paste0('Point pattern ', nam, '\nCell', namIn))
            points(coords(shifted$ppp), col = Cols[marks(shifted$ppp, drop = FALSE)$gene], pch = ".", cex = Cex)
            counter <- counter + 1
        }
    }
    addLegend(Cols, Shift <- c(counter %% Ceils[1], counter %/% Ceils[1]))
    if (colourBorder) {
        borderCols <- unique(unlist(borderCols))
        names(borderCols) <- unVals
        addLegend(borderCols, Shift + c(1, 0), Pch = 5, Main = borderColVar, Cex = 0.7)
    }
    if (warnPosition) {
        text(Ceils[1] / 2, -0.4, labels = "Cells not on original location but sorted by expression!")
    }
    invisible()
}
#' @importFrom spatstat.geom affine.owin affine.ppp
toUnitSquare <- function(owin, ppp, Shift) {
    shrink <- 1 / rep(max(diff(owin$xrange), diff(owin$yrange)), 2)
    vec <- -c(owin$xrange[1], owin$yrange[1]) * shrink + Shift
    list(owin = affine.owin(owin, diag(shrink), vec), ppp = affine.ppp(ppp, diag(shrink), vec))
}
