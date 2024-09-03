#' Plot the n cells with highest abundance of a feature
#' @description After testing for within-cell patterns, it may be useful to
#' look at the cells with the most events for certain genes. These are plotted
#' here, but the spatial location of the cells in the point pattern is lost! 
#' The choice and ranking of cells is one of decreasing gene (pair) expression.
#'
#' @param obj A hyperframe, or an object containing one
#' @param features The features to be plotted, a character vector
#' @param nCells An integer, the number of cells to be plotted
#' @param Cex The point expansion factor
#' @param borderColVar The variable to colour borders of the cell
#' @param borderCols Colour palette for the borders
#' @param warnPosition A boolean, should a warning be printed on the
#'  image that cells are not in their original location?
#' @param summaryFun A function to summarize the gene-cell table in case multiple genes are plotted.
#' Choose "min" for cells with the highest minimum, or "sum" for highest total expression
#' of the combination of genes
#' @param plotNuclei A boolean, should nuclei be added?
#' @param nucCol A character string, the colour in which the nucleus' boundary is plotted
#' @param ... Additional arguments, currently ignored
#' @inheritParams plotExplore
#'
#' @return Plots cells with highest expression to the plotting window, returns invisible
#' @export
#' @importFrom spatstat.geom plot.owin subset.ppp
#' @importFrom graphics points text
#'
#' @examples
#' example(addCell, "smoppix")
#' plotCells(hypFrame2, "gene1")
#' plotCells(hypFrame2, "gene1", borderColVar = "condition", nCells = 10)
plotCells <- function(
    obj, features = getFeatures(obj)[seq_len(3)], nCells = 100,
    Cex = 1.5, borderColVar = NULL, borderCols = rev(palette()), Mar = c(
        0.5, 0.1,
        0.75, 0.1
    ), warnPosition = TRUE, summaryFun = "min", 
    plotNuclei = !is.null(getHypFrame(obj)$nuclei), nucCol = "lightblue", ...) {
    if (!is.hyperframe(obj)) {
        obj <- getHypFrame(obj)
    }
    if(is.null(obj$owins)){
      stop("No windows present in object! Add them using addCells() first.")
    }
    summaryFun <- match.fun(summaryFun)
    colourBorder <- !is.null(borderColVar)
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar = Mar)
    features <- unique(unlist(lapply(features, sund)))
    stopifnot(is.hyperframe(obj), length(nCells) == 1, all(features %in%
        getFeatures(obj)))
    Cols <- makeCols(features, obj)
    # Show cells with highest minimum of events of the gene (pair)
    tablesCell <- lapply(obj$ppp, function(p) {
        tab <- table(marks(p[marks(p, drop = FALSE)$gene %in% features, ], drop = FALSE)[, c("gene", "cell")])
        apply(tab, 2, summaryFun, simplify = FALSE)
    })
    nCells <- min(nCells - 1 - colourBorder, length(ul <- unlist(tablesCell)))
    nthcell <- sort(ul, decreasing = TRUE)[min(nCells, length(ul))]
    tablesCell <- lapply(tablesCell, function(x) {
        x[x >= nthcell & names(x) != "NA"]
    })
    # Randomly subset to meet require cell number, starting from cells with low expression
    if ((nCellsEffective <- sum(lls <- vapply(tablesCell, FUN.VALUE = double(1), length))) > nCells) {
        idCells <- keepCells <- rep(seq_along(tablesCell), times = lls)
        ul <- unlist(tablesCell)
        counter <- 0
        while ((tooMany <- sum(idCells != 0) - nCells) > 0) {
            idMin <- which(ul <= (nthcell - counter))
            if (length(idMin) >= tooMany) {
                idCells[sample(idMin, tooMany)] <- 0
            } else {
                idCells[idMin] <- 0
            }
            counter <- counter + 1
        }
        tablesCell <- lapply(seq_along(tablesCell), function(x) {
            tablesCell[[x]][idCells[keepCells == x] == x]
        })
    } else {
        nCells <- nCellsEffective
    }
    names(tablesCell) <- rownames(obj)
    # Readjust to number of cells really used
    if (colourBorder) {
        if (borderColVar %in% getEventVars(obj)) {
            unVals <- unique(unlist(lapply(obj$ppp, function(x) unique(marks(x, drop = FALSE)[[borderColVar]]))))
            names(borderCols)[seq_along(unVals)] <- unVals
            borderCols <- lapply(seq_along(tablesCell), function(x) {
                borderCols[marks(obj$ppp[[x]])[[borderColVar]][match(
                    names(tablesCell[[x]]),
                    marks(obj$ppp[[x]])$cell
                )]]
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
    Ceils <- splitWindow(nCells + 1 + colourBorder)
    plot(
        type = "n", x = c(0, Ceils[1]), y = c(0 - warnPosition * 0.12, Ceils[2]),
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = FALSE
    )
    for (i in seq_along(tablesCell)) {
        nam <- names(tablesCell)[i]
        ppp <- subset.ppp(obj[[nam, "ppp"]], gene %in% features)
        for (j in seq_along(tablesCell[[nam]])) {
            namIn <- names(tablesCell[[nam]])[j]
            shifted <- toUnitSquare(obj[[nam, "owins"]][[namIn]], 
                                    ppp = subset.ppp(ppp, cell == namIn), 
                                    Shift = shiftVec(counter, Ceils[1]),
                                    nuclei = if(plotNuclei) obj[[nam, "nuclei"]][namIn])
            plot.owin(shifted$owin, add = TRUE, border = borderCols[[i]][[j]])
            if(plotNuclei){
              for(nn in shifted$nuclei){
                plot.owin(nn, add = TRUE, border = nucCol)
              }
            }
            # text(x = Shift[1], y = Shift[2], labels = paste0('Point pattern
            # ', nam, '\nCell', namIn))
            points(coords(shifted$ppp),
                col = Cols[marks(shifted$ppp, drop = FALSE)$gene],
                pch = ".", cex = Cex
            )
            counter <- counter + 1
        }
    }
    addLegend(Cols, shiftVec(counter, Ceils[1]))
    if (colourBorder) {
        borderCols <- unique(unlist(borderCols))
        names(borderCols) <- borderCols #unVals[[seq_along(borderCols)]]
        addLegend(borderCols, shiftVec(counter + 1, Ceils[1]),
            Pch = 5, Main = borderColVar,
            Cex = 0.7
        )
    }
    if (warnPosition) {
        text(Ceils[1] / 2, -0.12, labels = "Cells not in original location but sorted by expression!", cex = 0.75)
    }
    invisible()
}
#' @importFrom spatstat.geom affine.owin affine.ppp
toUnitSquare <- function(owin, ppp, Shift, nuclei) {
    shrink <- 1 / rep(max(diff(owin$xrange), diff(owin$yrange)), 2)
    vec <- -c(owin$xrange[1], owin$yrange[1]) * shrink + Shift
    list(owin = affine.owin(owin, diag(shrink), vec), 
         ppp = affine.ppp(ppp, diag(shrink),vec),
         nuclei = if(!is.null(nuclei)) {
           lapply(nuclei, affine.owin, diag(shrink), vec)
         })
}
shiftVec <- function(counter, Ceil) {
    c(counter %% Ceil, counter %/% Ceil)
}
